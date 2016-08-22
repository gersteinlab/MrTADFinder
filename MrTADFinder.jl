using HDF5; 
using JLD;
using MAT;
using Graphs;
using DataFrames;
#using Winston
using CurveFit;
using Distributions;

###################################################################################

function read_WG_contact_map(input_file,N);

	table=readdlm(input_file);
	A=sparse(table[:,1],table[:,2],table[:,3],N,N);
	tmp=A-spdiagm(diag(A));
	A=A+tmp';
	out_name="./all_contacts.jld";
	save(out_name,"interaction",A);

	return A;
end

function extract_chr(A,chr2bins,chr_num);
	st=1+chr2bins[1,chr_num];
	ed=1+chr2bins[2,chr_num];
	A_chr=A[st:ed,st:ed];
	return A_chr;
end

function get_expect_vs_d_v2(contacts,chr2bins,Nall);

	f_dall=zeros(Nall);
	tt_dall=zeros(Nall);

	intra_sum=zeros(24);

	contacts[isnan(contacts)]=0;

	S=sum(contacts);

	Neff=sum(sum(contacts,1).>0);

	Neff_intra=zeros(24);

	for chr_num=1:24

		display(chr_num);
	
		W=extract_chr(contacts,chr2bins,chr_num);
		W=full(W);
		W[isnan(W)]=0;
		dark_bins=find(sum(W,1).==0);
		N=size(W,1);
		f_d=zeros(N);
		tt_d=zeros(N);

		Neff_intra[chr_num]=N-length(dark_bins);
		intra_sum[chr_num]=sum(W);

		for d=0:N-1

			cd=diag(W,d);
			x=collect(1:N-d);
			y=collect(1+d:N);
			#d =1 means 1 vs 2, 2 vs 3.....
			is_okx=zeros(size(x));
			is_oky=zeros(size(y));
			for k=1:length(x)
				is_okx[k]=x[k] in dark_bins;
				is_oky[k]=y[k] in dark_bins;			
			end
			iz=find((1-is_okx).*(1-is_oky).>0);
			f_d[d+1]=sum(cd[iz]);
			tt_d[d+1]=length(iz);
		end

		f_dall[1:length(f_d)]=f_dall[1:length(f_d)]+f_d;
		tt_dall[1:length(f_d)]=tt_dall[1:length(f_d)]+tt_d;
	

	end


	expect_d=f_dall./tt_dall;
	#if we cf. this aggregated one with just chromosome 1, the number for d=a few are very consistent..
	inter_sum=S-sum(intra_sum);

	inter_chr_expect=inter_sum./(Neff^2-sum(Neff_intra.^2));


	return expect_d,inter_chr_expect;

end

function get_expect_vs_d(W);

	W[isnan(W)]=0;

	dark_bins=find(sum(W,1).==0);

	N=size(W,1);
	
	f_d=zeros(N);
	tt_d=zeros(N);

	#use the one with dark bins...
	#sum of all f_d are the reads..	
	#tt_d are the number of loci pairs separated by a distance d
	#f_d are the number of contacts between loci pairs separated by a distance d
	#f_d./tt_d is expectation..
	for d=0:N-1
		#display(d);
		cd=diag(W,d);
		x=collect(1:N-d);
		y=collect(1+d:N);
		#d =1 means 1 vs 2, 2 vs 3.....
		is_okx=zeros(size(x));
		is_oky=zeros(size(y));
		for k=1:length(x)
			is_okx[k]=x[k] in dark_bins;
			is_oky[k]=y[k] in dark_bins;			
		end
		iz=find((1-is_okx).*(1-is_oky).>0);
		f_d[d+1]=sum(cd[iz]);
		tt_d[d+1]=length(iz);
	end

	expect_d=f_d./tt_d;
	expect_d[isnan(expect_d)]=0;
	#NB. we average out the coverage for all loci in this calculation..

	return expect_d;

end

#a new fct for fitting expect_d. a general fct fitting everything may not work...
#y=Kx^-gamma;
function fit_expect_d(expect_d);

	x=collect(1:length(expect_d));

	iz=find(expect_d.==0)[1];
	x1=x[1:iz-1];
	y1=expect_d[1:iz-1];

	#x1=x[expect_d.>0];
	#y1=expect_d[expect_d.>0];
	#y[y.==0]=eps();
	#fit=curve_fit(PowerFit,x,y);
	tmp=linear_fit(log10(x1),log10(y1));#this line is the same as powerfit	

	#we cannot do that..too large contribution from the leading points.
	#tmp=linear_fit(log10(x_bin),log10(y_bin));
	
	gamma=tmp[2];
	K=10^tmp[1];
	return gamma, K;

end

function get_f_d_by_fitting(W,expect_d);

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));#what's f_W? it's a generation of ones(size(W));

	x=collect(1:length(expect_d));
	gamma,K=fit_expect_d(expect_d);
	expect_d2=K*x.^gamma;
	
	#this step is added for Ren's data..
	#expect_d2[1]=expect_d[1];
	#miss this thought..

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=expect_d2[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	f_W=f_W/sum(f_W)*N_eff.^2;

	return f_W;

end

function get_null_polymer(W,f_W,err_threshold);

	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);

	coverage_real=sum(W,2);
	invC=1./coverage_real;
	invC=invC[:,1];
	invC[isinf(invC)]=0;
	invC=diagm(invC);


	coverage_est=sum(W,2)/sqrt(sum(W));#sum(W) is 2N
	coverage_est[dark_bins]=0;

	#y=1./coverage_est;
	iy=coverage_est;
	
	tmp=f_W*iy;
	y_new=invC*tmp;
	coverage_est_new=1./y_new;
	coverage_est_new[isinf(coverage_est_new)]=0;
	#coverage_est_new=coverage_est_new/sum(coverage_est_new)*sqrt(sum(W)); this is not the right normalization..
	nm=sum((coverage_est_new*coverage_est_new').*f_W);
	coverage_est_new=coverage_est_new/sqrt(nm)*sqrt(sum(W));

	err=sum(abs(coverage_est_new-coverage_est));

	while err>err_threshold;
		display(err);
		coverage_est=coverage_est_new+0;
		iy=coverage_est;
		tmp=f_W*iy;
		y_new=invC*tmp;
		coverage_est_new=1./y_new;
		coverage_est_new[isinf(coverage_est_new)]=0;
		nm=sum((coverage_est_new*coverage_est_new').*f_W);
		coverage_est_new=coverage_est_new/sqrt(nm)*sqrt(sum(W));
		#coverage_est_new=coverage_est_new/sum(coverage_est_new)*sqrt(sum(W));
		err=sum(abs(coverage_est_new-coverage_est));
	end

	E_W=(coverage_est_new*coverage_est_new').*f_W;

	return coverage_est_new,E_W;
end

###################################################################################
function get_optimal_partition(W,E_W,res);
	
	Q=W-E_W*res;
	Optimal=zeros(size(Q));
	#traceback_aux=Array(Any,size(Q));
	traceback_aux=zeros(size(Q));
	for i=1:size(Q,1);
		Optimal[i,i]=Q[i,i];
		#traceback_aux[i,i]=[1];
		traceback_aux[i,i]=1;
	end
	N=size(Q,1);
	x=collect(1:N);

	for L=2:N;
		display(L);
		st,ed=get_all_substrings(x,L);
		st_pts=collect(1:L);
		cut_pts=collect(1:L+1:L^2);
		ed_pts=collect(L^2-L+1:L^2);
		for k=1:length(st)
			ix=collect(st[k]:ed[k]);
			Z=repmat(ix,1,L)';
			#pairs=get_decomposition(ix);
			possibility=zeros(L);
			
			for p=1:L-1
				#possibility[p]=Optimal[pairs[p,1][1],pairs[p,1][end]]+Optimal[pairs[p,2][1],pairs[p,2][end]];
				possibility[p]=Optimal[Z[st_pts[p]],Z[cut_pts[p]]]+Optimal[Z[cut_pts[p]+L],Z[ed_pts[p]]];
			end
			#possibility[1:L-1]=Optimal[sub2ind((N,N),Z[st_pts[1:end-1]],Z[cut_pts[1:end-1]])]+Optimal[sub2ind((N,N),Z[cut_pts[1:end-1]+L],Z[ed_pts[1:end-1]])];

			possibility[L]=sum(Q[ix,ix]);
			Optimal[st[k],ed[k]]=maximum(possibility);
			#traceback_aux[st[k],ed[k]]=find(possibility.==maximum(possibility));
			traceback_aux[st[k],ed[k]]=indmax(possibility);
		end
	end

	active_partition=[];
	final_partition=[];
	push!(active_partition,collect(1:N));

	while ~isempty(active_partition)
		partition=active_partition[1];
		partition_st=partition[1];
		partition_ed=partition[end];
		partition_len=partition_ed-partition_st+1;
		#pt=traceback_aux[partition_st,partition_ed][1];
		pt=traceback_aux[partition_st,partition_ed];
		if pt.==partition_len
			push!(final_partition,partition);
			shift!(active_partition);
		else
			push!(active_partition,partition[1:pt]);
			push!(active_partition,partition[pt+1:end]);
			shift!(active_partition);
		end
	end

	final_assign=zeros(Int,N);
	for i=1:length(final_partition)
		x=sort(final_partition[i]);
		final_assign[x[1]:x[end]]=i;
	end

	return final_assign;

end

function get_all_substrings(x,L);
	st=collect(1:x[end]-L+1);
	ed=st+L-1;
	return st,ed;
end

function get_decomposition(ix);
	m=length(ix)-1;
	pairs=cell(m,2);
	for pt=1:m
		pairs[pt,1]=ix[1:pt];
		pairs[pt,2]=ix[pt+1:end];
	end
	return pairs;
end



###################################################################################
#un-assigned bins are darks bin, not playing role...
function optimize_TADs_modlouvain(W,E_W,res,order=1);
#B is the generalized modularity matrix, and sW is sum of all elts in W, ~ to 2M
#perform the step of louvain
#order=1, forward
#order=0, random
	
	N=size(W,1);
    i_no_dark=find(sum(abs(W),1).>0);
    N_no_dark=length(i_no_dark);

    B=W-E_W*res;
    Bcompact=B[i_no_dark,i_no_dark];
    Wcompact=W[i_no_dark,i_no_dark];
    sW=sum(W);

    #(assign,Q,Brenorm)=iterate_TADs_modlouvain(Bcompact,sW,order);
    (assign,Q,Brenorm)=iterate_TADs_modlouvain_v2(Bcompact,sW,order);
    transfer=sparse(collect(1:size(assign,1)),assign,ones(size(assign)));

    keep_doing=1;

    while keep_doing==1
        #(tmp_assign,tmp_Q,tmp_Brenorm)=iterate_TADs_modlouvain(Brenorm,sW,order);
        (tmp_assign,tmp_Q,tmp_Brenorm)=iterate_TADs_modlouvain_v2(Brenorm,sW,order);
        tmp_transfer=sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        if isequal(tmp_transfer,speye(length(tmp_assign)))
        	keep_doing=0;
        	#Brenorm, Wrenorm are optimal.
        else
	        transfer=transfer*sparse(collect(1:size(tmp_assign,1)),tmp_assign,ones(size(tmp_assign)));
        	Brenorm=tmp_Brenorm+0;
        	Q=tmp_Q+0;
        end
    end

    (u,v)=findn(transfer);
    iu=sortperm(u);
    tmp_assign=v[iu];
    final_assign=zeros(Int,N);
    final_assign[i_no_dark]=tmp_assign;

    #pick up the zeros....dark bins...where should they belong?
    #look at two sides, if both sides are the same, merge them, if not, keep dark...
    #if the beginning is dark, i.e. loc=1. we kick it out
    #if the end is dark, i.e. loc[end:end+span-1]=0, we kick it out...
    (loc,span)=get_chunks_v2(final_assign,1);#
    #to break things into chunks, the number of chunks > number of domains above..
    #because 0 can be inserted into 2 sides with the same domains..
    if final_assign[1].==0
    	loc=loc[2:end];
    	span=span[2:end];
    end

    if final_assign[loc[end]].==0
    	loc=loc[1:end-1];
    	span=span[1:end-1];
    end

    for i=1:length(loc)
    	if final_assign[loc[i]].==0
    		if final_assign[loc[i]-1]==final_assign[loc[i]+span[i]];
    			final_assign[loc[i]:loc[i]+span[i]-1]=final_assign[loc[i]-1];
    		end
    	end
    end

    return final_assign, Q;

end

#This is a modiifed version of Louvain, update is done by comparing with nearest neighbors
#if order is 0, update is in random order..testing for a single iteration, it seems order is better than random
#if we remove dark bins before, the nearest neigthbor will be bins across the dark bins..
#if we don't remove dark bins, multiple dark bins in the middle will be hard to merge..
#I think in general, by removing the dark bins, the resultant Q is a bit higher. since dark bins have no contibution
#to Q at all theoretically, the result means the heuristic work better without dark bins, 
#for the optimal choice, removing dark bins wouldn't change the results, but make the matrice a bit smaller..

function iterate_TADs_modlouvain_v2(Bcompact,sW,order);

    Nb=size(Bcompact,1);

	if Nb>1
		sigma=collect(1:Nb);
		u=collect(1:Nb);
		if order==0
			u=u[randperm(Nb)];
		end

    	gain=1;Niter=1;

	    while (gain==1)
    	    gain = 0;
        	for j=1:Nb
        		#display(j);
            	x=u[j];
            	spin=sigma[x];
            	if x==1
            		spin_f=sigma[x+1];
            		spin_b=sigma[x];
            	elseif x==Nb;
            		spin_f=sigma[x];
            		spin_b=sigma[x-1];
            	else 
            		spin_f=sigma[x+1];
            		spin_b=sigma[x-1];
            	end
            	c=Bcompact[x,:];
           		c[x]=0;#this is important step to make sure the deltaQ is right 
            	neighbors_spin=sigma;
            	DeltaQ=-sum(c'.*(sigma.==spin))+full(sparse(neighbors_spin,[1 for dd=1:Nb],vec(c)));
            	#the 2nd term sum over the components from each community in advance
            	#1st term, the effect of getting rid of the original spin contribution..
            	#note the dim of DeltaQ is the number of communities
            	spin_choice=[spin_b spin spin_f];
            	id=indmax(DeltaQ[spin_choice]);#choose 1 out of 3...
            	id=spin_choice[id];
            	new_spin=id;
            	if (new_spin!=spin)&(DeltaQ[id].>0);
                	gain=1;
                	sigma[x]=new_spin;
            	end
        	end
        	Q=compute_modularity(sigma,Bcompact,sW);
        	@printf("iteration %d - sum = %f %d communities \n",Niter,Q,length(unique(sigma)))
        	Niter = Niter + 1
    	end

	    Q=compute_modularity(sigma,Bcompact,sW);

	    sigma2=relabel_communities(sigma);
    	usigma=sort(unique(sigma));
    	N_renorm=length(usigma);

  		Bcompact_renorm=zeros(N_renorm,N_renorm);

	    for i=1:N_renorm
	    	#display(i);
    	    for j=1:N_renorm
        	    Bcompact_renorm[i,j]=sum(sum(Bcompact[sigma.==usigma[i],sigma.==usigma[j]]));
        	end
    	end

	    Q2=compute_modularity(collect(1:N_renorm),Bcompact_renorm,sW);
    	@printf("step - Q = %f %d communities \n",Q2,length(unique(sigma)))
    	#println(Q2), verified to be the same as Q...
    else

		sigma2=[1];
		Q2=Bcompact/sW;
		Bcompact_renorm=Bcompact;

	end

    return sigma2, Q2, Bcompact_renorm;

end

function relabel_communities(sigma)
    u=unique(sigma);
    u=sort(u);
    sigma_new=zeros(size(sigma));
    for i=1:length(u)
        iz=findin(sigma,u[i]);
        sigma_new[iz]=i;
    end
    sigma_new=round(Int64,sigma_new);
    return sigma_new;
end


function compute_modularity(sigma,Brenorm,sW);
    COMu = unique(sigma);
    Q=0;
    for k=1:length(COMu)
        id = find(sigma.==COMu[k]);
        Q=Q+sum(Brenorm[id,id]);
    end
    Q=Q/sW;
    return Q;
end


#id is the starting loc of a chunk, and d is the length it spans..
function get_chunks_v2(a,singleton=0);
	# adopt from a matlab code by Jiro Doke;
	 a                 = [NaN; a; NaN];
	 b                 = diff(a);
	 b1                = b;  # to be used in fullList (below)
	 ii                = trues(size(b));
	 ii[b.==0] = false;
	 b[ii]             = 1;
	 c                 = diff(b);
	 id                = find(c.==-1);

	 #Get single-element chunks also
	 if singleton.==1
	 	b1[id]          = 0;
	 	ii2             = find(b1[1:end-1]);
	 	d               = vcat(find(c.==1) - id + 1, ones(length(ii2)));
	 	id              = [id;ii2];
	 	v=sortperm(id);
	 	id=sort(id);
	 	#(id,tmp)        = sort(id);
	 	d               = d[v];
	 else 
	 	d               = find(c.==1) - id + 1;
	 end

	 return id,d;
end

#this following code is the basic TAD calling code, if we adopt a convention to save the output at a particular 
#folder, we can use the folder to do a few downstream work
function get_high_confidence_boundaries_and_domains(W,E_W,res,num_trial,sig_cut);

	Z=zeros(size(W));
	all_bdd_rec=zeros(Int,size(W,1)+1,num_trial);
	all_final_assign=zeros(Int,size(W,1),num_trial);
	for x=1:num_trial;
		display(x);
		final_assign_x, Q1=optimize_TADs_modlouvain(W,E_W,res,0);
		all_final_assign[:,x]=final_assign_x;
		for cc=1:maximum(final_assign_x);
			idcc=find(final_assign_x.==cc)[1];
			all_bdd_rec[idcc,x]=1;
		end
		if final_assign_x[end].>0
			all_bdd_rec[end,x]=1;
		end
	end
	bdd_prob_score=mean(all_bdd_rec,2);
	#the actual modularity of the consensus domain is in fact lower that one the domains in one-trial
	#but the boundaries are more confident...
	
	tmp=find(bdd_prob_score.>=sig_cut);
	consensus_bdd=zeros(Int,size(bdd_prob_score));
	consensus_bdd[tmp]=1;
	consensus_domains=cumsum(consensus_bdd)[1:end-1];
	unassign_score=mean(all_final_assign.==0,2)
	i_unassign=find(unassign_score.>=sig_cut);
	consensus_domains[i_unassign]=0;
	TADs=consensus_domains+0;
	(a,b)=hist(TADs,collect(-.5:maximum(TADs)+1));
	b=b[2:end];
	#make sure all domains have pos number of reads..it common that some bins with zero
	#entry in diagonal form a single domains...
	aux=zeros(size(b));
	for j=1:length(b);
		iz=find(TADs.==j);
		if sum(W[iz,iz]).==0
			TADs[iz]=0;
			aux[j]=1;
		end
	end

	TADs_final=relabel_communities(TADs)-1;

	bdd_aux=sign(TADs_final);
	bdd_aux=[bdd_aux;0]+[0;bdd_aux];

	bdd_prob_score=bdd_prob_score.*sign(bdd_aux);

    return bdd_prob_score,bdd_aux,TADs_final;

    #bdd_prob_score.*(bdd_aux==2) to get the 2-sided boundaries
    #bdd_prob_score.*(bdd_aux==1) to get the 1-sided boundaries

end

#################################################################################################
#this code just report, no more filtering..
function report_domains(chr2bins,bin2loc,chr_num,TADs_final)

	u,v=get_chunks_v2(TADs_final,1);
	TAD_st=Int64[];
	TAD_ed=Int64[];

	for i=1:length(u);
		if TADs_final[u[i]]>0
			push!(TAD_st,u[i])
			push!(TAD_ed,u[i]+v[i]-1);
		end
	end
    
    TADs_list=DataFrame(chr=ASCIIString[],domain_st=Int64[],domain_ed=Int64[],domain_st_bin=Int64[],domain_ed_bin=Int64[],idx=Int64[]);

    st=chr2bins[:,chr_num][1]+1;#the bin count from zero in the stored file.
    ed=chr2bins[:,chr_num][2]+1;#we here shift it..
    chr2all_bin=collect(st:ed);
    #this will be the array mapped by elements in the chr of interest.

    if chr_num<=22
        chr_string=string("chr",string(chr_num));
    elseif chr_num==23
        chr_string=string("chr","X");
    elseif chr_num==24 
        chr_string=string("chr","Y");
    end;

    TAD_st_bin=chr2all_bin[TAD_st];
    TAD_ed_bin=chr2all_bin[TAD_ed];   

    for i=1:size(TAD_st_bin,1)
        loc1=bin2loc[2,TAD_st_bin[i]];
        loc2=bin2loc[3,TAD_ed_bin[i]];
        push!(TADs_list,[chr_string loc1 loc2 TAD_st_bin[i] TAD_ed_bin[i] i]);
    end

    return TADs_list
end

#we assume a folder TADs_loc that store the jld output of all chr
function concatenate_TADs_diff_chr(TADs_loc,TADs_file_prefix);
	max_TAD=0;
	all_TADs=[];
	for chr_num=1:24
		display(chr_num);
		TADs_file=TADs_loc*TADs_file_prefix*"_"*change_chr(chr_num)*".jld";
		d=load(TADs_file);
		TADs=d["TADs"];
		TADs[TADs.>0]=TADs[TADs.>0]+max_TAD;
		all_TADs=[all_TADs;TADs];
		max_TAD=maximum(TADs);
	end
	all_TADs=[all_TADs;0];

	return all_TADs;
	#this is the genome-wide bin to TAD map..
end


#function report_boundaries(consensus_domains,chr_num,bin2loc);

#	xxx1=zeros(Int,length(consensus_domains)+1);
#	for cc=1:maximum(consensus_domains);
#		idcc=find(consensus_domains.==cc)[1];
#		xxx1[idcc]=1;
#	end

#	xxx2=zeros(Int,length(consensus_domains)+1);
#	for cc=1:maximum(consensus_domains);
#		idcc2=find(consensus_domains.==cc)[end]+1;
#		xxx2[idcc2]=1;
#	end

#	xxx=sign(xxx1+xxx2);
    
#   chr_string=change_chr(chr_num);
#    TADs_boundaries=DataFrame(chr=ASCIIString[],Bdd=Int64[]);

#    i_bdd=find(xxx.>0);
#    i_chr=find(bin2loc[1,:].==chr_num-1);
#    bin_st=bin2loc[2,i_chr];
#    bin_ed=bin2loc[3,i_chr];
    
#    for i=1:length(i_bdd);
#    	push!(TADs_boundaries,[chr_string bin_st[i_bdd[i]]]);
#    end
#    if consensus_domains[end].>0
#		push!(TADs_boundaries,[chr_string bin_ed[end]]);
#	end

#    return TADs_boundaries;

#end

function generate_TADs_bed(TAD_list,out_file);

	X=[TAD_list[:chr] TAD_list[:domain_st] TAD_list[:domain_ed]];
	writedlm(out_file,X);

end

function show_mat(bins2modules);
    
    L=length(bins2modules);
    iz=find(bins2modules.>0);
    s2=bins2modules[iz];
   
    Z=(broadcast(-,s2,s2').==0)+0;
    
    allTADs=extend_mat(Z,iz,L);
    #imagesc(Z)

    return allTADs;

end

#the bed file has only chr num, st pos and end pos..
function read_TADs_bed(input_file,bin2loc,chr2bins,chr_num)
	all_Ren_TADs=readtable(input_file,header=false,separator='\t');
	#info=matread("/gpfs/scratch/fas/gerstein/ky26/Hi-C_processed/hg18_bin_info.mat");
	#chr2bins=int(info["chr2bins"]);
	#bin2loc=info["bin2loc"];
	chr=change_chr(chr_num);
	i_ren=find(all_Ren_TADs[:1].==chr);
	Ren_TADs_chr=all_Ren_TADs[i_ren,:];
	if chr_num.>1
		tmp=chr2bins[2,chr_num-1]+1;#tmp+1 correct the problem in previous round..
	else 
		tmp=1;
	end

	Ren_TADs_chr[:x4]=tmp+round(Int64,floor(Ren_TADs_chr[:x2]/40000)+1);
	Ren_TADs_chr[:x5]=tmp+round(Int64,floor(Ren_TADs_chr[:x3]/40000)+1);
	Ren_TADs_chr[:x6]=collect(1:size(Ren_TADs_chr,1));
	rename!(Ren_TADs_chr,:x1,:chr)
	rename!(Ren_TADs_chr,:x4,:domain_st_bin)
	rename!(Ren_TADs_chr,:x5,:domain_ed_bin)
	rename!(Ren_TADs_chr,:x6,:idx);
	rename!(Ren_TADs_chr,:x2,:domain_st)
	rename!(Ren_TADs_chr,:x3,:domain_ed)
	return Ren_TADs_chr;
end

#for unassigned bins, like dark ones, bins2modules will be zeros
function TADs_list_to_bins(TADs_list,chr2bins);
    
    chr_num=change_chr(TADs_list[:chr][1]);
    stc=chr2bins[:,chr_num][1]+1;#the bin count from zero in the stored file.
    edc=chr2bins[:,chr_num][2]+1;#we here shift it..
    bins2modules=zeros(Int,edc-stc+1,1);
    size_chr=length(bins2modules);

    chr2all_bin=collect(stc:edc);
    #this will be the array mapped by elements in the chr of interest.

    for i=1:size(TADs_list,1);
        stm=find(chr2all_bin.==TADs_list[i,4])[1];
        edm=find(chr2all_bin.==TADs_list[i,5])[1];
        bins2modules[stm:edm]=TADs_list[i,6];
    end

    return bins2modules;

end

function get_bdd_loc(is_bdd,chr_num,bin2loc)

	bin_size=bin2loc[3,1]-bin2loc[2,1]+1;
	i_bdd=find(is_bdd);
	i_pick=find(bin2loc[1,:].==chr_num-1);
	st=bin2loc[2,i_pick];
	ed=bin2loc[3,i_pick];
	bdd_loc=[st ed[end]];
	bdd_loc_array=zeros(Int,bdd_loc[end]);
	L=ed[end];
	bdd_loc=bdd_loc[i_bdd];
	for j=1:length(bdd_loc);
		iz=collect(bdd_loc[j]-Int(bin_size/2):1:bdd_loc[j]+Int(bin_size/2));
		iz=iz[(iz.>0).*(iz.<=L)];
		bdd_loc_array[iz]=1;
	end
	
	return bdd_loc,bdd_loc_array;
end

#use this code for boundary, which is defined by 2 bins...
#the first and last bdd therefor should not fit..
function get_local_interaction_strength(bdd,W,win_size);
	ib=collect(bdd[1]-1:bdd[1]);
	ia=collect(bdd[2]:bdd[2]+1);
	ia=ia[ia.>win_size-1];
	ib=ib[ib.<size(W,1)-win_size]+1;
	int_strength=(sum(W[ib,ib])+sum(W[ia,ia]))/(sum(W[ib,ia])+sum(W[ia,ib]));
	return int_strength;
end

function MI_two_partitions(a1,a2);

	iz=find(a1+a2.>0);
	a1=a1[iz];
	a2=a2[iz];
    m1=maximum(a1);
    m2=maximum(a2);
    (u1,v1)=hist(a1,-.5:m1+.5)
    (u2,v2)=hist(a2,-.5:m2+.5)
    (u12,v12,c12)=hist2d([a1 a2],-.5:m1+.5,-.5:m2+.5);
    P1=v1/sum(v1)
    P2=v2/sum(v2)
    P12=c12/sum(c12)
    
    H1=P1.*log2(P1);
    H1[isnan(H1)]=0;
    H1=-sum(H1)

    H2=P2.*log2(P2);
    H2[isnan(H2)]=0;
    H2=-sum(H2)

    Z=P12./(P1*P2');
    S=P12.*log2(Z);
    S[isnan(S)]=0;

    MI=sum(S);
    MI_norm=2*MI/(H1+H2);

    return MI, MI_norm

end

function swap_TADs(b2m);
	u,v=get_chunks_v2(b2m,1);
	r=randperm(length(u));
	b2m_r=zeros(size(b2m));
	count_pos=1;
	count_d=1;
	for i=1:length(u)
		j=r[i];
		if b2m[u[j]]>0
			b2m_r[count_pos:count_pos+v[j]-1]=count_d;
			count_d=count_d+1;
		else
			b2m_r[count_pos:count_pos+v[j]-1]=0;
		end
		count_pos=count_pos+v[j];
	end
	return b2m_r;
end

function swap_boundaries(is_bdd);
	is_bdd_r=BitArray(size(is_bdd));
	is_bdd_r[:]=false;
	u=find(is_bdd.>0);
	v=diff(u);
	v_r=v[randperm(length(v))];
	is_bdd_r[1]=true;
	count_pos=1;
	for i=1:length(v_r);
		count_pos=count_pos+v_r[i];
		is_bdd_r[count_pos]=true;
	end
	return is_bdd_r;
end


function extend_mat(Z,iz,L);
    (u,v)=ind2sub(size(Z),find(Z.!=0));
    w=Z[find(Z)];
    #w=nonzeros(Z);
    u=iz[u];
    v=iz[v];
    Z_extend=sparse(u,v,w,L,L);
    Z_extend=full(Z_extend);
    return Z_extend;
end

#it look at the enrichment of contacts over a polymer model.
#0 means no enrichment. value are -log10(P) based on a Poisson distribution.
function get_P_value_observed_vs_expected(W,E_polymer);
	all_P=ones(size(W));
	iz=find(sum(W,2).>0);
	for i=1:length(iz);
		display(i);
		for j=i:length(iz);
			if E_polymer[iz[i],iz[j]].>0 && W[iz[i],iz[j]].>E_polymer[iz[i],iz[j]]
				dd=Poisson(E_polymer[iz[i],iz[j]]);
				lw=cdf(dd,floor(W[iz[i],iz[j]]));
				uw=cdf(dd,ceil(W[iz[i],iz[j]]));
				delta=W[iz[i],iz[j]]-floor(W[iz[i],iz[j]]);
				cc=1-(lw+(uw-lw)*delta);
				all_P[iz[i],iz[j]]=cc;
				all_P[iz[j],iz[i]]=cc;
			end
		end
	end
	max_log_P=-floor(log10(minimum(all_P[all_P.>0])))+1;
	enrich=-log10(all_P);
	enrich[isinf(enrich)]=max_log_P;
	
	return enrich;
end

function get_null_d_local(W);

	N=size(W,1);
	E_W_local=zeros(size(W));
	for i=1:round(Int,ceil(N/2));
		z=W[i,:];
		zf=flipdim(z,2);
		z1=[zf[1:N-2*i+1]' z];
		#z1=[zeros(1,N-2*int(i)+1) z];
		#z2=[fliplr(z) zeros(1,N-2*int(i)+1)];
		z2=[zf z[2*i:end]'];
		#z1[N-i+1]=z2[N-i+1];
		z_null=(z1+z2)/2;
		st=N-2*i+2;
		ed=2*N-2*i+1;
		z_null=z_null[st:ed];
		E_W_local[i,:]=z_null;
	end
	for i=round(Int,ceil(N/2))+1:N;
		z=W[i,:];
		zf=flipdim(z,2);
		#z1=[z zeros(1,2*int(i)-N-1)];
		z1=[z flipdim(z[1:2*i-N-1]',2)];
		#z2=[zeros(1,2*int(i)-N-1) fliplr(z)];
		z2=[z[1:2*i-N-1]' flipdim(z,2)];
		#z1[i]=z2[i];
		z_null=(z1+z2)/2;
		st=1;
		ed=N;
		z_null=z_null[st:ed];
		E_W_local[i,:]=z_null;
	end

	return E_W_local;
	#off set by trace(W)/2

end


function get_f_d(W,expect_d)

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=expect_d[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	
	f_W=f_W/sum(f_W)*N_eff^2;
	return f_W;

end

