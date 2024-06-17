#! /home/xxia/usr/bin/python3

import os, sys, itertools
# from markers_4_unique_parents import Main
# from fileProcessing import Marker

def arg(argv):
    #comondline = '%s -sh %s -v %s -i %s -c %s --address %s -w %s -o %s &'%(pyscriptpath, rscriptpath, VCFfile, targetgeno, return_dict['checkout'], return_dict['email'], savedir, outfile)

    import argparse
    parser = argparse.ArgumentParser(description="snpEff piplines", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-sh', '--shscript', help="the different shell scripts to process different piplines", type=str, required=True)
    parser.add_argument('-f', '--vcffile', help="the vcffile that upload by users", type=str, required=True)
    parser.add_argument('-i', '--input', help="the species inputed by users", type=str, required=True)
    parser.add_argument('-c', '--checkout', help="the paramenters that corresponding to different output formats", type=str, nargs='*', default=None)
    parser.add_argument('--address', help="Email address specified if online tools are used.", type=str, default='FALSE')
    parser.add_argument('-w', '--savedir', help="the dirname that save the results", type=str, default=None)
    parser.add_argument('-o', '--output', help="Output information of final results", type=str, required=True)

    #upload genome need paramenters
    # parser.add_argument('-g', '--genome', help="the uploaded genome that be used to build database", type=str, default=None)
    # parser.add_argument('-e', '--gene', help="the uploaded gff/gtf/genebank file that be used to build database", type=str, default=None)

    args = parser.parse_args()
    # print(args)
    # print(args.shscript)
    #species, VCF, paramenters, savedir
    # para=""
    print(args.checkout)
    if args.checkout==[]:
        args.checkout="v"
    # para=' \-'.join(args.checkout)
    # print(para)
    # para2=args.checkout.replace(",", " \\-")

    main(args.shscript, args.input,  args.vcffile, args.checkout, args.savedir, args.output, args.address)

def main(script, species, vcf, paramemter, dirname, logfile, address):
    # print(args)
    # # sys.argv[3] == 'True':# sys.argv[3] == "False":# sys.argv[3] == "raise_Error":# sys.argv[3] == "raise_sys_Error"
    # comondline = 'nohup %s %s %s %s %s > %s 2>&1 &'%(script, species, vcf,paramemter,dirname,logfile)
    comondline = 'source %s %s %s %s %s > %s'%(script, species, vcf, paramemter, dirname, logfile)
    print(comondline)

    os.system(comondline)

    result_file=dirname+"/Result/Effresult.snp.eff.vcf"
    # print(result_file)
    send_file=dirname+"/Effresult.zip"

    if not os.path.exists(result_file):
        #file not exist, pipline not run
        #email: /home/lzhang/Django_SNPeff/project/static/SNPEFF/ptms/send_mail.py
        os.system("/home/xxia/usr/bin/python3 /home/lzhang/Django_project/send_mail_test/send_mail.py %s %s %s" % (address, send_file, 'raise_sys_Error'))
        sys.exit()

    elif not os.path.getsize(result_file):
        #pipline run but not get result 
        os.system("/home/xxia/usr/bin/python3 /home/lzhang/Django_project/send_mail_test/send_mail.py %s %s %s" % (address, send_file, 'raise_Error'))
        os.system("rm -rf %s/Result"%(dirname))
        sys.exit()
    else:
        os.system("/home/xxia/usr/bin/python3 /home/lzhang/Django_project/send_mail_test/send_mail.py %s %s %s" % (address, send_file, 'True'))
        # os.system("rm -rf %s/Result"%(dirname))
        sys.exit()


if __name__ == '__main__':
    arg(sys.argv)


# /home/lzhang/Django_SNPeff/project/static/SNPEFF/snpEff_file/Eff_now/run_snp.py -sh /home/lzhang/Django_project/send_mail_test/VCFfile.vcf -i grasscarp -c detail --address 2211696720@qq.com  -w /home/lzhang/Django_project/send_mail_test -o /home/lzhang/Django_project/send_mail_test/out