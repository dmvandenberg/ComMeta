from Bio import SeqIO
import re
import os


def check_path(path):
    """checks path"""
    try:
        os.remove(path)
    except:
        pass


def make_infoline(feature):
    """makes a string with info about a feature for the gff file"""
    if feature.type == "CDS":
        info_line = "ID={};locus_tag={};prodct={}".format(feature.qualifiers["protein_id"][0],
                                                          feature.qualifiers["locus_tag"][0],
                                                          feature.qualifiers["product"][0])
    elif feature.type == "misc_feature":
        info_line = "note={}".format(feature.qualifiers["note"][0])
    else:
        info_line = "locus_tag={};product={}".format(feature.qualifiers["locus_tag"][0],
                                                     feature.qualifiers["product"][0])

    return info_line


def to_gff(gbk, path):
    """extracts info and writes .gff file"""
    check_path(path)
    for record in SeqIO.parse(gbk, "genbank"):
        with open(path, "a") as gfffile:
            gfffile.write("##sequence-region {} 1 {}\n".format(record.id, len(record.seq)))
            for feature in record.features:
                if feature.type != "gene" and \
                        feature.type != "source" and \
                        feature.type != "assembly_gap":
                    if feature.type == "CDS":
                        try:
                            info_line = make_infoline(feature)
                            location = str(feature.location).split("](")[0][1:].split(":")
                            strand = str(feature.location).split("](")[1][0]
                            gfffile.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t0\t{}\n".format(record.id, feature.ref,
                                                                                      feature.type, location[0],
                                                                                      location[1], strand, info_line))
                        except KeyError:
                            pass
    with open(path, 'a') as gfffile:
        gfffile.write("##FASTA")
    return path


def to_faa(gbk, path):
    """extracts info and writes .faa file"""
    check_path(path)
    for record in SeqIO.parse(gbk, "genbank"):
        with open(path, "a") as faafile:
            for feature in record.features:
                if feature.type != "gene" and \
                        feature.type != "source" and \
                        feature.type != "assembly_gap":
                    if feature.type == "CDS":
                        try:
                            seq = feature.qualifiers["translation"][0]
                            seq = re.sub("(.{60})", "\\1\n", seq, 0, re.DOTALL)[:-1]
                            faafile.write(">{} {}\n{}\n".format(feature.qualifiers["protein_id"][0],
                                                                feature.qualifiers["product"][0],
                                                                seq))
                        except:
                            pass
    return path


def to_fna(gbk, path):
    """extracts info and writes .fasta file"""
    check_path(path)
    for record in SeqIO.parse(gbk, "genbank"):
        with open(path, "a") as fnafile:
            seq = re.sub("(.{70})", "\\1\n", str(record.seq), 0, re.DOTALL)
            fnafile.write(">{}\t{}\n{}\n".format(record.id,
                                                 record.description,
                                                 seq))
    return path
