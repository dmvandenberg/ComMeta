# necessary imports
import os
import console


def mummerplot(indir, subject, algorithm, outdir):
    """runs promer or nucmer and mummerplot to create a dotplot"""
    print "\tRunning {}:".format(algorithm)
    for f in os.listdir(indir):
        if f.endswith(".fasta") and f != subject:
            outname = outdir + "/" + subject.rsplit(".", 1)[0] + "+" + f.rsplit(".", 1)[0]
            os.system('{} {} {} --prefix {}'.format(algorithm, indir + '/' + subject, indir + '/' + f, outname))
            os.system("mummerplot {} --png --prefix {} --layout".format(outname + '.delta', outname))

            if outname.rsplit('/', 1)[1] + '.png' in os.listdir(outdir):
                os.system('rm {} && rm {} && rm {} && rm {}'.format(outname + '.delta', outname + '.rplot',
                                                                    outname + '.fplot', outname + '.gp'))
                print "\t\tCreated file {}".format(outname + '.png')


def run(indir, outdir, algorithm, subject):
    """main method for calling functions"""
    outdir = console.mk_dirs(outdir + '/dotplots/' + algorithm)
    if outdir[0]:
        mummerplot(indir, subject, algorithm, outdir[1])
