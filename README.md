# GPS is a clock problem

An interactive explainer. Your phone finds you by timing radio signals from satellites, but its clock is cheap. Change the clock error and watch the range circles stop meeting. Then let the receiver solve for its own clock and watch them snap back together.

**One thing you learn:** GPS needs one more satellite than the dimensions you are locating in (4 in real 3D, 3 in this 2D toy), because the extra satellite is spent on fixing your phone's clock, not on position.

## Run it

Open `index.html` in any current browser. There are no dependencies, no build step, no network calls, no sign-in and no API keys.

## Files

- `index.html` - the interactive page (canvas rendering, controls, explanations)
- `solver.js` - the mechanism, written from scratch: pseudorange model, Levenberg-Marquardt least squares for (x, y, clock bias), and HDOP from the inverse of J^T J

## What is computed live

Every circle radius, estimate, error figure, clock estimate and HDOP is computed in the browser from the current inputs (receiver position, satellite positions, clock error, timing noise, satellite count). Nothing is prerecorded.

## Host it (public link)

- GitHub Pages: push these files to a public repo, then Settings > Pages > deploy from `main` / root.
- Or drag the folder onto Netlify Drop / Cloudflare Pages.

## Honest simplifications

The map is a toy scale (about 3,000 km wide, satellites drawn close). Real GPS satellites orbit around 20,000 km up and the solve is in 3D. The mechanism is the same.
