
<!-- PROJECT SHIELDS -->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]



<!-- PROJECT LOGO -->
  <h3 align="center">The ICA Toolkit</h3>

  <p align="center">
    An easy and bioinformatic oriented package to find the independent component most related to a given factor!
    <br />
    <a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ"><strong> Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/j-solor/ICA_Toolkit/blob/main/02_TCGA/Example_TCGA.Rmd">View Demo</a>
    ·
    <a href="https://github.com/j-solor/ICA_Toolkit/issues">Report Bug</a>
    ·
    <a href="https://github.com/j-solor/ICA_Toolkit/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
        </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

The ICA toolkit aims to give the ordinay bioinformatician the option to use Independent Component analysis (ICA) with ease.

You can use this package in many ways:
* Find the most robust number of components to do ICA
* Find the most correlated component to a given fator.
* Either way, This package provides the tools to also explore all the components together, or particular components in depth for association.

At the moment the best source of information of this package is the [Demo R notebook](https://github.com/j-solor/ICA_Toolkit/blob/main/02_TCGA/Example_TCGA.Rmd) using TCGA data.

As everyone like plots their own way, the function here provided return ggplot object you can easily modify and combine  to suit your needs

The method here provided is based in the analysis performed in [R Nicolle et al 2020](https://doi.org/10.1016/j.ebiom.2020.102858)

For now the only ICA algorithm implemented inside the ICA_Toolkit is [JADE](https://www.jstatsoft.org/article/view/v028i06), used through the [JADE R package](https://cran.r-project.org/web/packages/JADE/index.html) 
<p align="right">(<a href="#top">back to top</a>)</p>

### Built With

This section should list any major frameworks/libraries used to bootstrap your project. Leave any add-ons/plugins for the acknowledgements section. Here are a few examples.

* R version 4.1.2 
* [tidyverse](https://www.tidyverse.org/)
* [JADE](https://cran.r-project.org/web/packages/JADE/index.html)
* [Hmisc](https://cran.r-project.org/web/packages/Hmisc/index.html)
* [pheatmap](https://cran.r-project.org/web/packages/pheatmap/)
* [scales](https://scales.r-lib.org/)
* [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
* [ggplotify](https://cran.r-project.org/web/packages/ggplotify/index.html)
* [patchwork](https://patchwork.data-imaginist.com/)

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Jacobo Solorzano  - jacobosolorzano96@gmail.com

Project Link: [https://github.com/j-solor/ICA_Toolkit](https://github.com/j-solor/ICA_Toolkit)

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- FUTURE IMPROVEMENTS -->
## Future improvements
* Fix Corr_by_Ct once Gemdecan from tpm is fixed 
* Give this a Package structure.
* Add bootstrapping to select the most robust number of components as an option.
* Add function to parse and get the best components.
* add other ICA algorithms·
* Little changes specified in each function documentation.
* Review documentation.


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/j-solor/ICA_Toolkit.svg?style=for-the-badge
[contributors-url]: https://github.com/j-solor/ICA_Toolkit/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/j-solor/ICA_Toolkit.svg?style=for-the-badge
[forks-url]: https://github.com/j-solor/ICA_Toolkit/network/members
[stars-shield]: https://img.shields.io/github/stars/j-solor/ICA_Toolkit?style=for-the-badge
[stars-url]: https://github.com/j-solor/ICA_Toolkit/stargazers
[issues-shield]: https://img.shields.io/github/issues/j-solor/ICA_Toolkit.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/j-solor/ICA_Toolkit/issues
[license-shield]: https://img.shields.io/github/license/j-solor/ICA_Toolkit?style=for-the-badge
[license-url]: https://github.com/j-solor/ICA_Toolkit/blob/master/LICENSE.txt
