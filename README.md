[![Licence](https://img.shields.io/github/license/Ileriayo/markdown-badges?style=for-the-badge)](./LICENSE)

# A Numerical Investigation of Stresses, Printing Efficiency, Printability, and Cell Viability in Nozzle Printheads for 3D Extrusion Bioprinting

## Introduction

3D (three-dimensional) extrusion bioprinting is a rapidly developing field in tissue engineering. Recent advances have ushered in a new stage of producing customized and bioengineered structures in regenerative medicine, pharmacokinetics, and basic cell biology studies. 3D extrusion bioprinters are the most widely used devices for various types of bioprinting devices because of their excellent cost-effectiveness and simple operation. They extrude bioinks (hydrogels that contain living cells) from needles to form filaments and scaffolds that can later be cross-linked and used to construct desired biostructures.

<p align="center">
  <img src="https://user-images.githubusercontent.com/111611023/203621835-3f76cf04-f0ea-444d-9015-1dc21c3deb09.png" alt="drawing" width="400"/>
</p>

During the extrusion, controlling stresses inside the needle is a significant factor in balancing printing resolution and cell viability. As experimentally observing stresses inside the syringe is complex, and testing thousands of bioink with different kinds of rheological properties is tedious and repetitive, the need to utilize numerical simulation to understand and optimize needle geometries, printing efficiency, and printability. Cell viability becomes an urgent task for 3D extrusion bioprinting devices. This research considers alginate-based bioinks since they are the most commonly used commercial bioinks due to low cost, biocompatibility, and facile gelation. Nonetheless, the numerical model can be easily modified to adapt to various kinds of bioinks.



<p align="center">
  <img src="https://user-images.githubusercontent.com/111611023/203621456-0d1ae19f-5506-49f7-9927-c2d86583a046.png" alt="drawing" width="500"/>
  <img src="https://user-images.githubusercontent.com/111611023/203620318-e4257e9f-9639-467e-94a6-4a78aad83c84.png" alt="drawing" width="350"/>
</p>

Current research focuses on utilizing OpenFOAM, an open-source software that implements computational fluid dynamics numerical solvers. Viscosity models, including power-law fluid and Herschel-Bulkley fluid, are considered. 2 types of needle shapes (cylindrical and tapered) and their corresponding shear stress distribution are investigated. The relationship among bioink's concentration, temperature, and apparent viscosity is evaluated. In addition, machine learning (supervised learning regression) is used to predict the cell viability inside the needle using both simulation and experimental results.
