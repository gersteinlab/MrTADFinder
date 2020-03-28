using DelimitedFiles
using LinearAlgebra;
using SparseArrays;

include("MrTADFinder.jl");

paras=[];
for x in ARGS;
	global paras
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
p=findfirst(isequal('='),res_str);
res=parse(Float64,res_str[p+1:end])

#chromosome 
pick_chr=paras[5];
pick_chr=parse(Int64,pick_chr);

#loc of output file
out_file=paras[6];
out_file2=out_file*".bddscore";

#Tunable parameters
num_trial=10;
sig_cut=.9;
err=1e-10;

#read binning information
println("reading binning information");
bins_info=readdlm(bins_file1);
chr2bins=bins_info[:,3:4];
chr2bins=chr2bins';
chr_num=bins_info[:,1];
chr_str=bins_info[:,2];
bin2loc=readdlm(bins_file2,Int64)';
N=chr2bins[2,pick_chr]+1;   # changed to chromosome specific 'N'
#N=chr2bins[2,end]+1;
#Nall=maximum(chr2bins[2,:]-chr2bins[1,:])+1;

bin_size=bin2loc[3,2]-bin2loc[3,1];

#load contact map
println("reading contact map");
contacts=read_generic_WG_contact_map(map_file,N);

W=extract_chr(contacts,chr2bins,pick_chr);
W2=W-diagm(0 => diag(W));
W=W2'+W;
W=Matrix(W);#full(W);
#use the next line to save the contact map in JLD format
#save("./all_contacts.jld","interaction",W);

#use this part if you want expect based on a genome-wide expectation
#xs_all,expect_d=get_expect_vs_d_WG_v0(contact,chr2bins,bin_size);

#expect based on single chr
println("estimating the null model");
xs_all,expect_d=get_expect_vs_d_single_chr_v0(W,chr2bins,bin_size);

#null model
f_W=get_f_W(W,expect_d);
c_est,E_W=get_null_polymer(W,f_W,err);

#TAD calling
println("start TAD calling");
bdd_prob_score,bdd_aux,TADs_final=get_high_confidence_boundaries_and_domains(W,E_W,res,num_trial,sig_cut);

#output
println("generate output");
TADs_list=report_domains(chr2bins,bin2loc,pick_chr,TADs_final);
generate_TADs_bed(TADs_list,out_file);

writedlm(out_file2,bdd_prob_score);