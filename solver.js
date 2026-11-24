/*
 * solver.js - the mechanism behind GPS positioning, written from scratch.
 *
 * Each satellite broadcasts "I am at s_i, and I sent this at time t".
 * The receiver subtracts that from its own clock's arrival time, giving a
 * "pseudorange":
 *
 *     rho_i = |p - s_i| + c*b          (distance + light-travel of the clock error)
 *
 * p is the unknown position, b is the unknown error of the receiver's clock.
 * Every satellite's circle is too big (or too small) by the SAME amount, c*b.
 * We solve for (x, y, b) by nonlinear least squares (Levenberg-Marquardt).
 * Distances here are in world pixels; the UI converts to km and microseconds.
 */
(function (root) {
  'use strict';

  var C_KM_PER_US = 0.299792458; // light travels ~0.3 km per microsecond

  // Solve A x = b by Gaussian elimination with partial pivoting. Returns null if singular.
  function solveLinear(A, b) {
    var n = b.length, M = [], i, r, c, k;
    for (i = 0; i < n; i++) M.push(A[i].concat([b[i]]));
    for (c = 0; c < n; c++) {
      var p = c;
      for (r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      if (Math.abs(M[p][c]) < 1e-9) return null;
      var tmp = M[c]; M[c] = M[p]; M[p] = tmp;
      for (r = c + 1; r < n; r++) {
        var f = M[r][c] / M[c][c];
        for (k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    var x = new Array(n);
    for (i = n - 1; i >= 0; i--) {
      var s = M[i][n];
      for (k = i + 1; k < n; k++) s -= M[i][k] * x[k];
      x[i] = s / M[i][i];
    }
    return x;
  }

  /*
   * sats: [{x,y}], rho: measured pseudoranges (same length)
   * estimateBias=false -> "trust my clock": unknowns x,y, b fixed at 0
   * estimateBias=true  -> real GPS: unknowns x,y,b
   */
  function solve(sats, rho, estimateBias) {
    var n = sats.length, m = estimateBias ? 3 : 2, i, j;
    if (n < m) return { ok: false, reason: 'underdetermined' };

    var x = 1.234, y = 0.567, b = 0;
    for (i = 0; i < n; i++) { x += sats[i].x / n; y += sats[i].y / n; }

    function cost(px, py, pb) {
      var c = 0;
      for (var q = 0; q < n; q++) {
        var f = Math.hypot(px - sats[q].x, py - sats[q].y) + pb - rho[q];
        c += f * f;
      }
      return c;
    }
    function build(px, py, pb) {
      var A = [], g = [];
      for (var a = 0; a < m; a++) { A.push(new Array(m).fill(0)); g.push(0); }
      for (var q = 0; q < n; q++) {
        var d = Math.max(Math.hypot(px - sats[q].x, py - sats[q].y), 1e-9);
        var row = [(px - sats[q].x) / d, (py - sats[q].y) / d];
        if (m === 3) row.push(1);
        var f = d + pb - rho[q];
        for (var a2 = 0; a2 < m; a2++) {
          g[a2] += row[a2] * f;
          for (var c2 = 0; c2 < m; c2++) A[a2][c2] += row[a2] * row[c2];
        }
      }
      return { A: A, g: g };
    }

    var lambda = 1e-2, c = cost(x, y, b), it, step = 1;
    for (it = 0; it < 150; it++) {
      var sys = build(x, y, b), improved = false;
      for (var tries = 0; tries < 14; tries++) {
        var Ad = sys.A.map(function (r, ri) {
          return r.map(function (v, ci) { return v + (ri === ci ? lambda * v + 1e-12 : 0); });
        });
        var d = solveLinear(Ad, sys.g.map(function (v) { return -v; }));
        if (d) {
          var nx = x + d[0], ny = y + d[1], nb = m === 3 ? b + d[2] : 0;
          var nc = cost(nx, ny, nb);
          if (nc < c) {
            x = nx; y = ny; b = nb; c = nc;
            lambda = Math.max(lambda / 3, 1e-9);
            step = Math.hypot.apply(null, d);
            improved = true;
            break;
          }
        }
        lambda *= 4;
      }
      if (!improved || step < 1e-10) break;
    }

    // Geometry quality: invert J^T J at the solution. If it can't be inverted, the
    // circles cannot pin the receiver down (satellites in a line, or circles that never meet).
    var fin = build(x, y, b), inv = [];
    for (j = 0; j < m; j++) {
      var e = new Array(m).fill(0); e[j] = 1;
      var col = solveLinear(fin.A, e);
      if (!col) return { ok: false, reason: 'singular', rms: Math.sqrt(c / n) };
      inv.push(col);
    }
    return {
      ok: true, x: x, y: y, b: b, iters: it,
      rms: Math.sqrt(c / n),
      hdop: Math.sqrt(inv[0][0] + inv[1][1])
    };
  }

  var api = { solve: solve, C_KM_PER_US: C_KM_PER_US };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.GPSSolver = api;
})(typeof window !== 'undefined' ? window : globalThis);
