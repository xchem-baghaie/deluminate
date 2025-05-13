import pandas as pd
from PIL import Image, ImageOps, ImageDraw, ImageFont

"""
Comparing virtual libraries to synthesized libraries.
"""


def image_grid(imgs, rows, cols):
    #assert len(imgs) == rows * cols

    w, h = imgs[0].size
    grid = Image.new("RGB", size=(cols * w, rows * h))

    for i, img in enumerate(imgs):
        grid.paste(img, box=(i % cols * w, i // cols * h))
    return grid


def generate_grid(
    in_filename: str = "/mnt/p/discovery chemistry/library synthesis/enumeration/library enumeration/Postsynthesis_Enumeration_Tracking.xlsx",
    fingerprint: str = "ECFP6",
    dimensionality_reduction: str = "UMAP",
    syn_syn_dir: str = "/mnt/p/Discovery Chemistry/People/Ryan/Tech Dev - Cheminformatics/Virtual Libraries/Chemical Space Comparisons/Synthesized vs Synthesized/",
    syn_vir_dir: str = "/mnt/p/Discovery Chemistry/People/Ryan/Tech Dev - Cheminformatics/Virtual Libraries/Chemical Space Comparisons/Synthesized vs Virtual/",
):
    df = pd.read_excel(in_filename, sheet_name="Postsynthesis Tracking")
    df = df.loc[
        (df["Ran Postsynthesis Enumeration"] == "Y")
        & (df["Ran Virtual Enumeration"] == "Y")
    ].reset_index(drop=True)
    df = df.sort_values(by=["Deck", "Library"]).reset_index()
    assert len(df) == 96

    lib_deck_dict = {df["Library"][i]: df["Deck"][i] for i in range(len(df))}
    deck_color_dict = {"DELcore": "#98D700", "DELflex": "#00A3EB"}

    syn_syn_img_ls = []
    syn_vir_img_ls = []
    for lib, deck in list(lib_deck_dict.items()):
        print(lib)
        ss = Image.open(
            f"{syn_syn_dir}{lib}/Static Plots/{dimensionality_reduction}/{fingerprint}/{lib}_Library Type.png"
        )
        ss = ImageOps.expand(ss, border=(25, 25, 25, 25), fill=deck_color_dict[deck])
        ss = ImageOps.expand(ss, border=(10, 10, 10, 10), fill="#FFFFFF")
        syn_syn_img_ls.append(ss)
        sv = ImageOps.expand(sv, border=(25, 25, 25, 25), fill=deck_color_dict[deck])
        sv = ImageOps.expand(sv, border=(10, 10, 10, 10), fill="#FFFFFF")
        syn_vir_img_ls.append(sv)
    syn_syn_grid = image_grid(syn_syn_img_ls, 6, 13)
    syn_vir_grid = image_grid(syn_vir_img_ls, 6, 13)
    syn_syn_grid.save(
        f"{syn_syn_dir}Synthesized_vs_Synthesized_{dimensionality_reduction}_{fingerprint}_per_library.png"
    )
    syn_vir_grid.save(
        f"{syn_vir_dir}Synthesized_vs_Virtual_{dimensionality_reduction}_{fingerprint}_per_library.png"
    )

    return None


def generate_grid_xbbs_only(
    in_filename: str = "/mnt/p/discovery chemistry/People/Ryan/Tech Dev - Cheminformatics/Virtual Libraries/Chemical Space Comparisons/Library Groups.csv",
    fingerprint: str = "ECFP6",
    dimensionality_reduction: str = "UMAP",
    syn_vir_dir: str = "/mnt/p/Discovery Chemistry/People/Ryan/Tech Dev - Cheminformatics/Virtual Libraries/Chemical Space Comparisons/XBBs Only/",
):
    df = pd.read_csv(in_filename).fillna('None')

    lib_deck_dict = {df["Library Group"][i]: df["Prioritization"][i] for i in range(len(df))}
    deck_color_dict = {"Productivity": "#00A3EB", "Size": "#98D700", "Both": "#F04178", "None":"#FFFFFF"}

    syn_vir_img_ls = []
    for lib, deck in list(lib_deck_dict.items()):
        print(lib)
        sv = Image.open(
            f"{syn_vir_dir}{lib}/Static Plots/{dimensionality_reduction}/{fingerprint}/{lib}_Library Type.png"
        )
        sv = ImageOps.expand(sv, border=(0,250,0,0), fill="#FFFFFF")
        sv = ImageOps.expand(sv, border=(25, 25, 25, 25), fill=deck_color_dict[deck])
        sv = ImageOps.expand(sv, border=(10, 10, 10, 10), fill="#FFFFFF")
        syn_vir_img_ls.append(sv)

    syn_vir_grid = image_grid(syn_vir_img_ls, 6, 13)

    syn_vir_grid.save(
        f"{syn_vir_dir}Synthesized_vs_Virtual_{dimensionality_reduction}_{fingerprint}_per_library.png"
    )

    return None

if __name__ == "__main__":
    generate_grid()
