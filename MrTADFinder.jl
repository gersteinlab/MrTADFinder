using HDF5; 
using JLD;
using MAT;
using Graphs;
using DataFrames;
#using Winston
using CurveFit;
using Distributions;

#text files in form of 3 columns..
function read_HiCPro_outputs(input_file_loc,bin_size);
	file_loc=input_file_loc*"/"*bin_size;
	filename=file_loc*"/reads_"*bin_size*"_iced.matrix";
	table=readtable(filename,separator='\t',header=false);
	N=maximum(table[:x1]);
	A=sparse(table[:x1],table[:x2],table[:x3],N,N);
	tmp=A-spdiagm(diag(A));
	A=A+tmp';
	#A=full(A);
	out_name=file_loc*"/contact.jld";
	save(out_name,"interaction",A);

	return A;

end

#NB. chr num starts with 1 to 25
function extract_chr(A,chr2bins,chr_num);
	st=1+chr2bins[1,chr_num];
	ed=1+chr2bins[2,chr_num];
	A_chr=A[st:ed,st:ed];
	return A_chr;
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
	y1=expect_d[1:iz-1]

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

function get_optimal_partition(W,E_W,gamma)
	
	Q=W-E*gamma;
	Optimal=zeros(size(Q));
	traceback_aux=Array(Any,size(Q));
	for i=1:size(Q,1);
		Optimal[i,i]=Q[i,i];
		traceback_aux[i,i]=[1];
	end
	N=size(Q,1);
	x=collect(1:N);

	for L=2:N;
		st,ed=get_all_substrings(x,L);
		for k=1:length(st)
			ix=collect(st[k]:ed[k]);
			pairs=get_decomposition(ix);
			possibility=zeros(L);
			for p=1:L-1
				possibility[p]=Optimal[pairs[p,1][1],pairs[p,1][end]]+Optimal[pairs[p,2][1],pairs[p,2][end]];
			end
			possibility[L]=sum(Q[ix,ix]);
			Optimal[st[k],ed[k]]=maximum(possibility);
			traceback_aux[st[k],ed[k]]=find(possibility.==maximum(possibility));
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
		pt=traceback_aux[partition_st,partition_ed][1];
		if pt.==partition_len
			push!(final_partition,partition);
			shift!(active_partition);
		else
			push!(active_partition,partition[1:pt]);
			push!(active_partition,partition[pt+1:end]);
			shift!(active_partition);
		end
	end
	return final_partition;

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




