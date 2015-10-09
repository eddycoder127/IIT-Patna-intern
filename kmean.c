#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define sqr(x) ((x)*(x))
#define MAX_CLUSTERS 20
#define MAX_ITERATIONS 400
#define BIG_double (INFINITY)

void fail(char *str)
{
	printf("%s\n",str);
	exit(-1);
}
double calc_distance(int dim, double *p1, double *p2)
{
	double distance_sq_sum=0;
	int ii;
	for(ii = 0;ii<dim;ii++)
		distance_sq_sum+=sqr(p1[ii]-p2[ii]);
	return distance_sq_sum;
}
void calc_all_distances(int dim,int n,int k,double *X,double *centroid, double *distance_output)
{
	int ii,jj;
	for(ii=0;ii<n;ii++)
		for(jj=0;jj<k;jj++)
		{
			distance_output[ii*k+jj]=calc_distance(dim,&X[ii*dim],&centroid[jj*dim]);
		}
}
double calc_total_distance(int dim,int n,int k,double *X,double *centroids,int *cluster_assignment_index)
{
	double tot_D=0;
	int ii;
	for(ii=0;ii<n;ii++)
	{
		int active_cluster=cluster_assignment_index[ii];
		if(active_cluster!=-1)
			tot_D+=calc_distance(dim,&X[ii*dim],&centroids[active_cluster*dim]);
	}
	return tot_D;
}
void choose_all_clusters_from_distances(int dim,int n,int k,double *distance_array,int *cluster_assignment_index)
{
	int ii,jj;
	for(ii=0;ii<n;ii++)
	{
		int best_index=-1;
		double closest_distance=BIG_double;
		for(jj=0;jj<k;jj++)
		{
			double cur_distance=distance_array[ii*k+jj];
			if(cur_distance<closest_distance)
			{
				best_index=jj;
				closest_distance=cur_distance;
			}
		}
		cluster_assignment_index[ii]=best_index;
	}
}
void calc_cluster_centroids(int dim,int n,int k,double *X,int *cluster_assignment_index,double *new_cluster_centroid)
{
	int cluster_member_count[MAX_CLUSTERS],ii,jj;
	for(ii=0;ii<k;ii++)
	{
		cluster_member_count[ii]=0;
		for(jj=0;jj<dim;jj++)
			new_cluster_centroid[ii*dim+jj]=0;
	}
	for(ii=0;ii<n;ii++)
	{
		int active_cluster=cluster_assignment_index[ii];
		cluster_member_count[active_cluster]++;
		for(jj=0;jj<dim;jj++)
			new_cluster_centroid[active_cluster*dim +jj]+=X[ii*dim+jj];
	}
	for(ii=0;ii<k;ii++)
	{
		if(cluster_member_count[ii]==0)
			printf("WARNING: Empty cluster %d \n",ii);
		for(jj=0;jj<dim;jj++)
			new_cluster_centroid[ii*dim+jj]/=cluster_member_count[ii];
	}
}
void get_cluster_member_count(int n,int k,int *cluster_assignment_index,int *cluster_member_count)
{
	int ii,jj;
	for(ii=0;ii<k;ii++)
		cluster_member_count[ii]=0;
	for(ii=0;ii<n;ii++)
		cluster_member_count[cluster_assignment_index[ii]]++;
}
void update_delta_score_table(int dim,int n,int k,double *X,int *cluster_assignment_cur,double *cluster_centroid,int *cluster_member_count,double *point_move_score_table,int cc)
{
	int ii,kk;
	for(ii=0;ii<n;ii++)
	{
		double dist_sum=0;
		for(kk=0;kk<dim;kk++)
		{
			double axis_dist=X[ii*dim+kk]-cluster_centroid[cc*dim+kk];
			dist_sum+=sqr(axis_dist);
		}
		double mult = ((double)cluster_member_count[cc]/(cluster_member_count[cc]+((cluster_assignment_cur[ii]==cc) ? -1 : +1)));
		point_move_score_table[ii*dim+cc]=dist_sum*mult;
	}
}
void perform_move(int dim,int n,int k,double *X,int *cluster_assignment,double *cluster_centroid,int *cluster_member_count,int move_point,int move_target_cluster)
{
	int cluster_old=cluster_assignment[move_point];
	int cluster_new=move_target_cluster;
	cluster_assignment[move_point]=cluster_new;
	cluster_member_count[cluster_old]--;
	cluster_member_count[cluster_new]++;
	if(cluster_member_count[cluster_old]<=1)
		printf("WARNING: Can't handle single-member clusters! \n");
	int ii;
	for(ii=0;ii<dim;ii++)
	{
		cluster_centroid[cluster_old*dim+ii]-=(X[move_point*dim+ii]-cluster_centroid[cluster_old*dim+ii])/cluster_member_count[cluster_old];
		cluster_centroid[cluster_new*dim+ii]+=(X[move_point*dim+ii]-cluster_centroid[cluster_new*dim+ii])/cluster_member_count[cluster_new];
	}
}
void cluster_diag(int dim,int n,int k,double *X,int *cluster_assignment_index,double *cluster_centroid)
{
	int cluster_member_count[MAX_CLUSTERS];
	get_cluster_member_count(n,k,cluster_assignment_index,cluster_member_count);
	printf("Final clusters \n");
	int ii;
	for(ii=0;ii<k;ii++)
		printf("cluster %d:    members: %8d,centroid(%.1f %.1f) \n",ii,cluster_member_count[ii],cluster_centroid[ii*dim+0],cluster_centroid[ii*dim+1]);
}
void copy_assignment_array(int n,int *src,int *tgt)
{
	int ii;
	for(ii=0;ii<n;ii++)
		tgt[ii]=src[ii];
}
int assignment_change_count(int n,int a[],int b[])
{
	int change_count=0;
	int ii;
	for(ii=0;ii<n;ii++)
		if(a[ii]!=b[ii])
			change_count++;
	return change_count;
}
void kmeans(int dim,double *X,int n,int k,double *cluster_centroid,int *cluster_assignment_final)
{
	double *dist=(double *)malloc(sizeof(double)*n*k);
	int *cluster_assignment_cur=(int *)malloc(sizeof(int)*n);
	int *cluster_assignment_prev=(int *)malloc(sizeof(int)*n);
	double *point_move_score=(double *)malloc(sizeof(double)*n*k);
	if(!dist||!cluster_assignment_cur||!cluster_assignment_prev||!point_move_score)
		fail("Error allocating dist arrays");
	calc_all_distances(dim,n,k,X,cluster_centroid,dist);
	choose_all_clusters_from_distances(dim,n,k,dist,cluster_assignment_cur);
	copy_assignment_array(n,cluster_assignment_cur,cluster_assignment_prev);
	double prev_totD=BIG_double;
	int batch_iteration=0;
	while(batch_iteration<MAX_ITERATIONS)
	{
		calc_cluster_centroids(dim,n,k,X,cluster_assignment_cur,cluster_centroid);
		double totD=calc_total_distance(dim,n,k,X,cluster_centroid,cluster_assignment_cur);
		if(totD>prev_totD)
		{
			copy_assignment_array(n,cluster_assignment_cur,cluster_centroid);
			calc_cluster_centroids(dim,n,k,X,cluster_assignment_cur,cluster_centroid);
			printf("   negative progress made on this step - iteration completed (%.2f) \n",totD-prev_totD);
			break;
		}
		copy_assignment_array(n,cluster_assignment_cur,cluster_assignment_prev);
		calc_all_distances(dim,n,k,X,cluster_centroid,dist);
		choose_all_clusters_from_distances(dim,n,k,dist,cluster_assignment_cur);
		int change_count=assignment_change_count(n,cluster_assignment_cur,cluster_assignment_prev);
		printf("%3d   %u    %9d   %16.2f    %17.2f\n",batch_iteration,1,change_count,totD,totD-prev_totD);
	       	fflush(stdout);
		if(change_count==0)
		{
			printf("no change made on this step - iteration completed \n");
			break;
		}
		prev_totD=totD;
		batch_iteration++;
	}
	cluster_diag(dim,n,k,X,cluster_assignment_cur,cluster_centroid);
	copy_assignment_array(n,cluster_assignment_cur,cluster_assignment_final);
	
	free(dist);
	free(cluster_assignment_cur);
	free(cluster_assignment_prev);
	free(point_move_score);
}
int main()
{
	double X[100000],cluster_centroid[100000];
	int n,dim,k,i,c,cluster_assignment_final[100000];
	long int f;
	i=0;
	scanf("%d%d%d",&n,&dim,&k);
	for(i=0;i<n*dim;i++)
		scanf("%lf",&X[i]);
	
			
	for(i=0;i<k*dim;i++)
	{
		cluster_centroid[i]=X[i];
	}
	for(i=0;i<n;i++)
	{
		cluster_assignment_final[i]=0;
	}
	kmeans(dim,X,n,k,cluster_centroid,cluster_assignment_final);
	for(i=0;i<n;i++)
	{
		printf("%d ",cluster_assignment_final[i]);
	}
	printf("\n");
	return 0;
}
