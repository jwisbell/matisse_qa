# import img2pdf
from glob import glob
from sys import argv
import numpy as np
from pypdf import PdfWriter, PdfReader
import matplotlib.pyplot as plt
from datetime import datetime


def convert_pngs_to_pdf(png_list, output_pdf):
    with open(output_pdf, "wb") as f:
        f.write(img2pdf.convert(*png_list[::-1]))


def create_title_page(output_path, title, subtitle, tau0, qc1, qc2, qc3):
    fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 size in inches (210mm x 297mm)
    ax.axis("off")  # No axes

    # Title
    plt.text(0.5, 0.75, title, ha="center", va="center", fontsize=24, weight="bold")

    # Subtitle
    plt.text(0.5, 0.65, subtitle, ha="center", va="center", fontsize=18)

    plt.text(0.5, 0.55, f"{tau0}", ha="center", va="center", fontsize=14)

    # Date
    plt.text(
        0.5,
        0.50,
        f"Created: {datetime.now().strftime("%B %d, %Y at %I:%M %p")}",
        ha="center",
        va="center",
        fontsize=12,
    )

    plt.text(
        0.5,
        0.45,
        f"{qc1}",
        ha="center",
        va="center",
        fontsize=12,
    )
    plt.text(
        0.5,
        0.40,
        f"{qc2}",
        ha="center",
        va="center",
        fontsize=12,
    )
    plt.text(
        0.5,
        0.35,
        f"{qc3}",
        ha="center",
        va="center",
        fontsize=12,
    )
    # Save as PDF
    plt.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close()


def merge_pdfs(pdf_paths, output_pdf):
    # Create a merger object
    writer = PdfWriter()
    print("writing pdf")

    for path in pdf_paths:
        reader = PdfReader(path)
        for page in reader.pages:
            writer.add_page(page)

    with open(output_pdf, "wb") as f_out:
        writer.write(f_out)


if __name__ == "__main__":
    script, fdir = argv

    # Get the list of PNG file paths
    png_files = np.sort(glob(f"{fdir}/*.png"))
    # Convert the PNGs to PDF
    convert_pngs_to_pdf(png_files, "./output.pdf")
