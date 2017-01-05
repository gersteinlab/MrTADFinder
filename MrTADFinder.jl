using JLD;
using DataFrames;
using CurveFit;
using Distributions;
using Interpolations;

###################################################################################
function read_generic_WG_contact_map(map_file,N);
	
	map=readdlm(map_file);
	I=map[:,1];
	I=round(Int64,I);
	J=map[:,2];
	J=round(Int64,J);
	K=map[:,3];
	W=sparse(I,J,K,N,N);
	W=full(W);

	return W;
end

##redundant function in HiC_spector..
function extract_chr(A,chr2bins,chr_num);
	st=1+chr2bins[1,chr_num];
	ed=1+chr2bins[2,chr_num];
	A_chr=A[st:ed,st:ed];
	return A_chr;
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

#this fct is the same as the one in HiC-spector
function get_expect_vs_d_single_chr_v0(W,chr2bins,bin_size);

	W=full(W);
	W[isnan(W)]=0;
	
	N=size(W,1);
	

	(u,v,w)=findnz(triu(W));
	d=float(v-u);
	d2=float(d);
	d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
	d3=d2*bin_size;

	#model = loess(log10(d3),log10(w),span=0.01);
	#the loess fct is rather slow, and fail to work at some matrices (not sure why), we have replaced it by a simpler method

	x=log10(d3);
	y=log10(w);

	xs,ys_smooth=local_smoothing(x,y);

	xs_all=collect(0:1.0:size(W,1)-1);xs_all[1]=1/3;
	xs_all=xs_all*bin_size;
	xs_all_aux=log10(xs_all);

	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all_aux);
		ik=find(xs.==xs_all_aux[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end	

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth,Gridded(Linear()));
	#itp=interpolate(knots,ys_smooth[ys_all.>0], Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=10.^ys_all;

	return xs_all, expect;

end

#this fct is the same as the one in HiC-spector
function get_expect_vs_d_WG_v0(contact,chr2bins,bin_size);

	#to find distance dependence, we should NOT iced the chr one by one.
	#because in genome-wide scale dependance, we should keep the contacts in same base

	all_d2=Float64[];
	all_w=Float64[];
	Ltmp=zeros(23);
	for chr_num=1:23
	
		#display(chr_num);
		W=extract_chr(contact,chr2bins,chr_num);
		W=full(W);
		W[isnan(W)]=0;

		N=size(W,1);
		
		(u,v,w)=findnz(triu(W));
		
		d=float(v-u);
		d2=float(d);
		d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
		
		all_d2=[all_d2;d2];
		all_w=[all_w;w];
		Ltmp[chr_num]=size(W,1);
	
	end

	all_d3=all_d2*bin_size;

	x=log10(all_d3);
	y=log10(all_w);

	xs,ys_smooth=local_smoothing(x,y);

	xs_all=collect(0:1.0:maximum(Ltmp)-1);xs_all[1]=1/3;
	xs_all=xs_all*bin_size;
	xs_all_aux=log10(xs_all);

	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all_aux);
		ik=find(xs.==xs_all_aux[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end	

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=10.^ys_all;

	return xs_all, expect;

end

#this fct is the same as the one in HiC-spector
function get_f_W(W,ys);

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));#what's f_W? it's a generation of ones(size(W));

	x=collect(1:N);

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=ys[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	f_W=f_W/sum(f_W)*N_eff.^2;

	return f_W;

end

###################################################################################
#this following code is the basic TAD calling code
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
	#make sure all domains have pos number of reads..it is common that some bins with zero
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

	#bdd_prob_score.*(bdd_aux==2) to get the 2-sided boundaries
    #bdd_prob_score.*(bdd_aux==1) to get the 1-sided boundaries

    return bdd_prob_score,bdd_aux,TADs_final;

end

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
           		c[x]=0;#this is important step to make sure the deltaQ is right ;

            	neighbors_spin=sigma;
            	DeltaQ=-sum(c'.*(sigma.==spin))+full(sparse(neighbors_spin,[1 for dd=1:Nb],squeeze(c',2)));
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
#this code is used in HiC-spector
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

function generate_TADs_output(TADs_list,out_file);
	writetable(out_file, TADs_list);
end

function generate_TADs_bed(TADs_list,out_file);
	X=[TADs_list[:chr] TADs_list[:domain_st] TADs_list[:domain_ed]];
	writedlm(out_file,X);
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


##############################################################################
#a few useful codes for downstream analysis

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

###############################################################

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

