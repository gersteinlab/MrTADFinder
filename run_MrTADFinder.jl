include("/home/fas/gerstein/ky26/Github/MrTADFinder/MrTADFinder.jl");

paras=[];
for x in ARGS;
	paras=vcat(paras,x);
end

#Input: 
#contact map 
map_file=paras[1];

#binning information
bins_file1=paras[2];
bins_file2=paras[3];

#user defined resolution parameter
res_str=paras[4];
p=search(res_str,"=")[end];
res=parse(Float64,res_str[p+1:end])

#chromosome 
pick_chr=paras[5];
pick_chr=parse(Float64,pick_chr);

#loc of output file
out_file=paras[6];

#Tunable parameters
num_trial=10;
sig_cut=.9;
err=1e-10;

#read binning information
bins_info=readdlm(bins_file1);
chr2bins=bins_info[:,3:4];
chr2bins=chr2bins';
chr_num=bins_info[:,1];
chr_str=bins_info[:,2];
bin2loc=readdlm(bins_file2,Int64)';
N=chr2bins[2,end]+1;
Nall=maximum(chr2bins[2,:]-chr2bins[1,:])+1;

bin_size=bin2loc[3,2]-bin2loc[3,1];

#load contact map
contacts=read_generic_WG_contact_map(map_file,N);

W=extract_chr(contacts,chr2bins,pick_chr);
W=full(W);
#use the next line to save the contact map in JLD format
#save("./all_contacts.jld","interaction",W);

#
#use this part if you want expect based on a genome-wide expectation
#xs_all,expect_d=get_expect_vs_d_WG_v0(contact,chr2bins,bin_size);

#expect based on single chr
xs_all,expect_d=get_expect_vs_d_single_chr_v0(W,chr2bins,bin_size);

#null model
f_W=get_f_W(W,expect_d);
c_est,E_W=get_null_polymer(W,f_W,err);

#TAD calling
bdd_prob_score,bdd_aux,TADs_final=get_high_confidence_boundaries_and_domains(W,E_W,res,num_trial,sig_cut);

#output
TADs_list=report_domains(chr2bins,bin2loc,pick_chr,TADs_final);
generate_TADs_bed(TADs_list,out_file);

