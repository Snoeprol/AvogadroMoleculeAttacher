import tkinter as tk
from . import parser

def generate():
    print("Crystal Atom:\t\t\t{}\nMolecule Attachment Atom:\t{}\nCrystal Attachment Atom:\t{}\nReplacement Atom:\t\t{}"\
    .format(e1.get(), e2.get(), e3.get(), e4.get()))
    parser.main(e1.get(),e2.get(),e3.get(),e4.get())

if __name__ == "__main__":
    master = tk.Tk()
    master.title("Meerio Crystal Attacher")
    tk.Label(master, 
            text="Crystal Atom (e.g. : Si):").grid(row=0)
    tk.Label(master, 
            text="Molecule Attachment Atom (e.g. : P):").grid(row=1)
    tk.Label(master,
            text = "Crystal Attachment Atom (e.g. : C):").grid(row=2)
    tk.Label(master,
            text = "Replacement Atom for Attachment atom (e.g. : O):").grid(row=3)

    e1 = tk.Entry(master)
    e2 = tk.Entry(master)
    e3 = tk.Entry(master)
    e4 = tk.Entry(master)
    e1.grid(row=0, column=1)
    e2.grid(row=1, column=1)
    e3.grid(row=2, column=1)
    e4.grid(row=3, column=1)

    tk.Button(master, 
            text='Quit', 
            command=master.quit).grid(row=5, 
                                        column=0, 
                                        sticky=tk.W, 
                                        pady=4)
    tk.Button(master, 
            text='Generate', command=generate).grid(row=5, 
                                                        column=1, 
                                                        sticky=tk.W, 
                                                        pady=4)
    tk.mainloop()
