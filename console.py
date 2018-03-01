# necessary imports
import argparse
import os
import shutil
from datetime import datetime
import gbk_parse
import identity
import dotplot
import comparaVis
import singleVis
import draw_venn

# program start time
t0 = datetime.now()

#########################################################################
# Description:  ComMeta                                                 #
# Version:      1.0                                                     #
# Dependencies: Circos, USearch, CompareM, Prodigal, Mummer3, Blast+    #
#               OrthoANIu, VennDiagram(R package)                       #
# Author:       D. van den Berg                                         #
# TBD: Accept fasta files as input as well (console.py)                 #
#########################################################################


def arguments():
    """arguments able to be parsed into the script"""
    parser = argparse.ArgumentParser(description="Welcome to ComMate!")

    # full workflow
    parser.add_argument('-i', '--indir',
                        required=True,
                        help="path to input directory (.gbff files)")
    parser.add_argument('-o', '--outdir',
                        required=True,
                        help="path to output directory")
    parser.add_argument('-s', '--subject',
                        required=True,
                        help="file name of subject genome (.fna)")
    parser.add_argument('-t', '--threads',
                        default=8,
                        type=int,
                        help="number of threads (default = 8)")
    parser.add_argument('--window',
                        type=int,
                        default=500,
                        help="window gc_skew and gc percentage (default = 500)")
    parser.add_argument('--kmer',
                        type=int,
                        default=1000,
                        help="minimum size for overlapping region visualisation (default = 1000)")
    parser.add_argument('--usearch',
                        default='/usr/local/bioinfo/usearch/usearch9',
                        help='path to usearch9 (default = "/usr/local/bioinfo/usearch/usearch9"')
    parser.add_argument('--treshold',
                        default=1e-6,
                        help="treshold for overlapping region e-value (default = 1e-6)")
    parser.add_argument('--no_ident',
                        action='store_true',
                        help="skip identity calculation")
    parser.add_argument('--no_dp',
                        action='store_true',
                        help="skip dotplot creation")
    parser.add_argument('--no_venn',
                        action='store_true',
                        help="skip venn creation")
    parser.add_argument('--no_scirc',
                        action='store_true',
                        help="skip single genome circos plot")
    parser.add_argument('--no_ccirc',
                        action='store_true',
                        help="skip comparative genome circos plot")
    parser.set_defaults(func=full_wflow)

    return parser.parse_args()


def mk_dirs(path):
    """Creates directory"""
    try:
        if os.path.exists(path):
            os.system('rm -rf {}'.format(path))
        os.makedirs(path)
        return True, path
    except:
        return "[ERROR]: path '{}' invalid\n".format(path)


def check_paths(indir, outdir, subject, usearch):
    """Checks if paths exists, creates them if they don't"""
    if os.path.exists(indir) and os.path.exists(outdir) and \
            subject in os.listdir(indir) and os.path.exists(usearch):
        return True
    elif not os.path.exists(usearch):
        print "[ERROR]: usearch not found"
    elif os.path.exists(indir) and subject not in os.listdir(indir):
        print "[ERROR]: subject file not found"
    elif not os.path.exists(indir):
        print "[ERROR]: input directory not found"
    elif not os.path.exists(outdir):
        path = mk_dirs(outdir)
        if path[0]:
            print "[MESSAGE]: output directory '{}' created".format(path[1])
            return True
        else:
            print "[ERROR]: something went wrong"


def parse_gbff(indir, outdir):
    """Generates .faa .fna and .gff for reference files"""
    outdir = mk_dirs("{}/generated".format(outdir))
    if outdir[0]:
        print ("[CONVERTING GBFF FILES]")
        for f in os.listdir(indir):
            if f.endswith(".gbff"):
                out_file = f.rsplit(".", 1)[0]
                base_path = outdir[1] + "/" + out_file
                gbk_parse.to_faa("{}/{}".format(indir, f), base_path + ".faa")
                gbk_parse.to_gff("{}/{}".format(indir, f), base_path + ".gff")
                gbk_parse.to_fna("{}/{}".format(indir, f), base_path + ".fasta")
                print ("\t{} converted".format(f))


def select(generated, outdir, dataset):
    """Places 4 genomes closest to subject in ./selected_data"""
    outdir = mk_dirs("{}/selected_data".format(outdir))
    for f in os.listdir(generated):
        if f.rsplit(".", 1)[0] in dataset:
            shutil.copyfile(generated + "/" + f, outdir[1] + "/" + f)
    return outdir[1]


def check_called(prokka_dir, subject, file_dir):
    """Check if prokka already called the subject file"""
    subname = subject.rsplit(".", 1)[0]
    try:
        if subname + ".faa" in os.listdir(prokka_dir) and \
                subname + ".gff" in os.listdir(prokka_dir):
            shutil.copyfile(prokka_dir + "/" + subname + ".gff", file_dir + "/" + subname + ".gff")
            shutil.copyfile(prokka_dir + "/" + subname + ".faa", file_dir + "/" + subname + ".faa")
        else:
            return False
    except:
        return False


def call_genes(subject, outdir, threads):
    """Calls genes on subject genome using prokka"""
    tmp_out = outdir + "/tmp/"
    subname = subject.rsplit(".", 1)[0].rsplit("/", 1)[1]
    print ("[MESSAGE]: Calling genes for {}".format(subject))
    try:
        os.system("prokka {} -o {} --prefix {} --metagenome --quiet --cpus {}".format(subject, tmp_out, subname,
                                                                                      threads))
        print ("[MESSAGE]: Gene calling done")
    except:
        print ("[ERROR]: Prokka gene calling failed")

    for out_file in os.listdir(tmp_out):
        if out_file.endswith(".faa") or out_file.endswith(".gff"):
            shutil.move(tmp_out + out_file, "{}/{}".format(outdir, out_file))

    shutil.rmtree(tmp_out)


def full_wflow(args):
    """call all the scripts in the correct order"""
    """Calls commands for all possible steps"""

    file_dir = args.outdir + "/generated"
    prokka_dir = args.outdir + "/prokka"

    # parse gbff to fasta faa and gff
    parse_gbff(args.indir, args.outdir)

    # copy subject fasta to generated for ANI
    shutil.copyfile(args.indir + "/" + args.subject, file_dir + "/" + args.subject)

    # get 4 genomes closest related to subject (ANI)
    identity.run(file_dir, args.outdir, args.threads, args.usearch, "ani")
    ani_dest = "{}/identity/ani_output/ani_output.txt".format(args.outdir)
    dataset = identity.get_dataset(ani_dest, args.subject)

    # make dataset directory
    selected_dir = select(file_dir, args.outdir, dataset)

    # called genes on subject fasta if not already existent
    prokka = check_called(prokka_dir, args.subject, selected_dir)
    if not prokka:
        call_genes(args.indir + "/" + args.subject, prokka_dir, args.threads)

    # place all required files in same folder
    shutil.copyfile(args.indir + "/" + args.subject, selected_dir + "/" + args.subject)
    for f in os.listdir(prokka_dir):
        shutil.copyfile(prokka_dir + "/" + f, selected_dir + "/" + f)

    # run AAI
    identity.run(selected_dir, args.outdir, args.threads, args.usearch, "aai")

    # get RBH's
    t1 = datetime.now()
    print ("[CALCULATING RECIPROCAL BEST BLAST HITS]")
    print ("=" * 57)
    os.system("seblastian.py -i {} -x faa -o {} -t {}".format(selected_dir, args.outdir + "/RBH", args.threads))
    print ("===duration: {}".format(datetime.now() - t1) + "=" * 30)

    # dotplots for nucmer & promer
    if not args.no_dp:
        t1 = datetime.now()
        print ("[GENERATING DOT PLOTS]")
        print ("=" * 57)
        dotplot.run(selected_dir, args.outdir, "nucmer", args.subject)
        dotplot.run(selected_dir, args.outdir, "promer", args.subject)
        print ("===duration: {}".format(datetime.now() - t1) + "=" * 30)

    # create single genome circos plot
    if not args.no_scirc:
        t1 = datetime.now()
        print ("[GENERATING GENOME CIRCOS PLOTS]")
        print ("=" * 57)
        singleVis.run(selected_dir, args.outdir, args.window)
        print ("===duration: {}".format(datetime.now() - t1) + "=" * 30)

    # create comparative genome circos plot
    if not args.no_ccirc:
        t1 = datetime.now()
        print ("[GENERATING COMPARATIVE CIRCOS PLOTS]")
        print ("=" * 57)
        comparaVis.run(selected_dir, args.outdir, args.subject, args.kmer, args.treshold)
        print ("===duration: {}".format(datetime.now() - t1) + "=" * 30)

    # draw venn diagram
    if not args.no_venn:
        t1 = datetime.now()
        print ("[GENERATING VENN DIAGRAM]")
        print ("=" * 57)
        draw_venn.run(selected_dir, args.outdir)
        print ("===duration: {}".format(datetime.now() - t1) + "=" * 30)


if __name__ == '__main__':
    """main method of the application"""
    args = arguments()
    paths = check_paths(args.indir, args.outdir, args.subject, args.usearch)

    if paths:
        args.func(args)
        print ("=" * 57)
        print ("total duration: " + str(datetime.now() - t0))
    else:
        quit()
