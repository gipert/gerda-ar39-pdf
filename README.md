<img src=".github/gerda-logo.png" align="left"  height="80"/>

# gerda-ar39-pdf

**[under development]** Ar39 pdfs for the GERDA experiment.

### Usage

Download the latest `gerda-ar39-lookup.tar.xz` from the [releases
page](https://github.com/gipert/gerda-ar39-pdf/releases). Extract the archive
in the repository root folder. You'll just need to include
`include/gerda_ar39_pdf.hpp` in order to be able to use the software:

```cpp
#include "include/gerda_ar39_pdf.hpp"

int main() {

    // the pdf is a function of germanium channel, energy (keV),
    // detector full charge-collection depth (FCCD) (mm), fully
    // dead-layer fraction (DLF)
    gerda::ar39_pdf(23, 123.4, 1.67, 0.44);

    return 0;
}
```
