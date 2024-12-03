import img2pdf
from glob import glob
from sys import argv
import numpy as np


def convert_pngs_to_pdf(png_list, output_pdf):
    with open(output_pdf, "wb") as f:
        f.write(img2pdf.convert(*png_list[::-1]))


script, fdir = argv

# Get the list of PNG file paths
png_files = np.sort(glob(f"{fdir}/*.png"))
# Convert the PNGs to PDF
convert_pngs_to_pdf(png_files, "./output.pdf")
