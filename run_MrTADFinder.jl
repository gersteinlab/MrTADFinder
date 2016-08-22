include("/home/fas/gerstein/ky26/Github/MrTADFinder/MrTADFinder.jl");

paras=[];
for x in ARGS;
	paras=vcat(paras,x);
end

#need a file of whole contact map, as output of ICED
#/home/fas/gerstein/ky26/scratch/Hi-C_data/Ren/IMR90_combine_out/hic_results/matrix/reads/iced/40000/
map_file=paras[1];
bins_file=paras[2];
res=paras[3];
chr_num=paras[4];

num_trial=10;
sig_cut=.9;
err=1e-10;

bins_info=readdlm(bins_file);
chr2bins=bins_info[:,3:4];
chr2bins=chr2bins';
chr_num=bins_info[:,1];
chr_str=bins_info[:,2];

N=chr2bins[2,end]+1;
contacts=read_WG_contact_map(map_file,N);
W=extract_chr(contacts,chr2bins,chr_num);

Nall=maximum(chr2bins[2,:]-chr2bins[1,:])+1;
expect_d,inter_chr_expect=get_expect_vs_d_v2(contacts,chr2bins,Nall);

f_W=get_f_d_by_fitting(W,expect_d);
c_est,E_W=get_null_polymer(W,f_W,err);

bdd_prob_score,bdd_aux,TADs_final=get_high_confidence_boundaries_and_domains(W,E_W,res,num_trial,sig_cut);

