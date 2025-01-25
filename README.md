# A simple CPA and DFA in Julia

This is a simple Julia implementation for performing a power correlation analysis attack. The code follows the structure outlined below:

- **Plaintext:** A matrix of shape `(X, 16)`, where `X` is the number of plaintexts sent, and `16` represents the division of each plaintext into 16 bytes.
- **Traces:** A matrix of shape `(X, N)`, where `X` is the number of captured traces (equal to the number of plaintexts sent), and `N` represents the number of points per trace.

## Example & Reference
The attack images provided in this repository are based on an example using `Coastalwhite`'s traces and plaintext. You can find his repository [here](https://github.com/coastalwhite/intro-power-analysis).

## Prerequisites
- Julia installed on your system
- Required Julia packages:
  - `Statistics`
  - `NPZ`
  - `Plots`
  - `ProgressBars`
  - `Printf`

## Coming soon
A julia code for a simple DFA