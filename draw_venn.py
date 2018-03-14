# necessary imports
from jinja2 import Environment, FileSystemLoader, Template
import os
from collections import OrderedDict
import console


def render_template(data, org_list, venn_dir, alg):
    """render values into the quad_venn template"""
    env = Environment(loader=FileSystemLoader('{}/templates'.format(os.path.abspath(__file__).rsplit("/",1)[0])))
    color_l = ["dodgerblue", "orchid3", "darkorange1", "seagreen3", "goldenrod1"]
    colors = ""
    for color in color_l[:alg[0]]:
        colors += '"{}",'.format(color)
    with open("{}/venn_conf.R".format(venn_dir), 'w') as f:
        r_template = env.get_template("venn.R")
        f.write(r_template.render(alg=alg[1], data=data, orgs=org_list, colors=colors[:-1],
                                  outname="{}/venn.tiff".format(venn_dir)))


def get_org_list(outdir):
    """makes dictionary with organism name and a assigned number"""
    org_dict = OrderedDict()
    number = 0
    for f in os.listdir("{}/RBH/output".format(outdir)):
        if f.startswith("unique"):
            number += 1
            org_dict[number] = f[7:-4]
    return org_dict


def get_areas(file_dir, org_dict):
    """gets total amount of genes in .faa files for each organism"""
    data = ""
    for f in os.listdir(file_dir):
        if f.endswith(".faa"):
            lines = open("{}/{}".format(file_dir, f), 'r').readlines()
            size = 0
            area = "area"
            for line in lines:
                if line.startswith(">"):
                    size += 1
            for key in org_dict:
                if org_dict[key] in f:
                    area += str(key)
            data += "{} = {},\n".format(area, size)
    return data


def get_info_string(file_dir, outdir, org_dict):
    """creates a string to put into R script template containing the info of overlap between genomes"""
    data = get_areas(file_dir, org_dict)
    for f in os.listdir("{}/RBH/output".format(outdir)):
        size = len(open("{}/RBH/output/{}".format(outdir, f), 'r').readlines())
        if f.startswith("unique"):
            continue
        else:
            comparison = "n"
            for key in org_dict:
                if org_dict[key] in f:
                    comparison += str(key)
        data += "{} = {},\n".format(comparison, size)
    return data


def get_alg(org_dict):
    """decides if a quad or a quintuple venn should be created based on number of organisms"""
    num = len(org_dict)
    alg = ""
    if num == 4:
        alg = "quad"
    elif num == 5:
        alg = "quintuple"
    return num, alg


def run(file_dir, outdir):
    """main method to call different functions"""
    venn_out = console.mk_dirs("{}/venn".format(outdir))
    org_dict = get_org_list(outdir)
    data = get_info_string(file_dir, outdir, org_dict)
    org_list = ""
    for key in org_dict:
        org_list += '"{}",'.format(org_dict[key])
    alg = get_alg(org_dict)
    render_template(data, org_list[:-1], venn_out[1], alg)
    os.system("Rscript {}".format(venn_out[1] + "/venn_conf.R"))
