# RationalizeRoots
RationalizeRoots is a Mathematica and Maple package for the rationalization of square roots.

## Installation

For both, Mathematica and Maple, all of the necessary source code is contained in a single file that you can find in this GitHub repository:

`RationalizeRoots.m` for Mathematica and `RationalizeRoots.mpl` for Maple.

In Mathematica, you can load the package via

`Get["<file location>/RationalizeRoots.m"]`

In Maple, you can load it via

`read "<file location>/RationalizeRoots.mpl"`

## Usage

In Mathematica:

```
RationalizeRoot[Sqrt[1-x^2-y^2]]
{{x -> 2t[1] / (1+t[1]^2+t[2]^2), y -> -(1-t[1]^2+t[2]^2) / (1+t[1]^2+t[2]^2)}}
```

## Contributors

[Marco Besier](marcobesier.com), Pascal Wasser, Stefan Weinzierl

## Further Reading

You can learn more about the functionality of the package and its mathematical background in the original [paper](https://arxiv.org/pdf/1910.13251.pdf).

## License

RationalizeRoots is released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).
