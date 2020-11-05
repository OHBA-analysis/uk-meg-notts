function [  ] = compute_global_env_netmat( Sin )

% [ res ] = compute_global_env_netmat( S )
%
% S.Ds=Ds
% S.num_iters=1000;
% S.roinets_protocol='symmetric'
% S.prefix='sf_'

Ds=Sin.Ds;

try roinets_protocol=Sin.roinets_protocol; catch roinets_protocol='symmetric'; end;
try centre_freq=Sin.centre_freq; catch centre_freq=13; end;

%% compute embedded cov matrices
S=[];
S.concat = [];
S.concat.protocol=roinets_protocol;
S.concat.embed.do=1;
S.concat.embed.centre_freq=centre_freq;
S.concat.embed.rectify=false;
S.concat.whiten=1;
S.concat.normalisation='voxelwise';
S.concat.pcadim=-1;
S.netmat_method=@netmat_cov;

[ state_netmats_cov ] = hmm_full_global_cov( Ds, S );

