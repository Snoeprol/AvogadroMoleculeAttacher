-------------PURPOSE OF CODE-----------------------

This project is used to attach molecules to crystals, using avogadro .xyz files.

-------------REQUIRMENTS---------------------------

To run this project, one needs python3 and the module numpy.
------------HOW TO RUN-----------------------------

1. Design your crystal and molecule in avogadro and put attachment sites on both the crystal and molecule.
You can put as many attachment sites on the crystal, but the molecule only needs one.
The attachment sites have to be different from the other atoms.
2. In the MoleculeAttacher folder navigate to the input folder and put in your molecule in the molecule folder
and the crystal in the crystal folder.
3. In your command prompt, move to the molcule attacher folder:
$cd MoleculeAttacher
4. Run the GUI.py file by typing:
$python -m src.GUI.py
5. An interface will open where you can put in your chosen atoms. Don't use quotation marks.
6. Press generate. The generated crystal with the molecules will be located in the output folder.# AvogadroMoleculeAttacher
