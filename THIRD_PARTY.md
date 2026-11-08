# Third-party software

`third_party/frlib` contains the required MATLAB runtime subset of
`frankpermenter/frlib`, with its original MIT license and README. The source
workspace records upstream commit `347f5b6201238cdf8b2ea11d5a0370190961589e`;
the manifest identifies the exact distributed bytes. This is not a claim that
the upstream Git object was fetched during this packaging step.

MATLAB, Symbolic Math Toolbox, YALMIP and MOSEK are external dependencies and
are not redistributed. Only new numerical searches require YALMIP and MOSEK.
The exact checks use native MATLAB and Symbolic Math Toolbox.
