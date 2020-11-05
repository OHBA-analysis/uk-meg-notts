function [ state_netmats ] = hmm_state_switch_mt_type( state_netmats, type, full_type )

% [ state_netmats ] = hmm_state_switch_mt_type( state_netmats, type, full_type )

for ss = 1:length(state_netmats),
    sub_state_netmats=state_netmats{ss};
    for k = 1:length(sub_state_netmats.state),
        data=[];

        S=[];
        S.netmats=sub_state_netmats.state{k};
        S.type=type;
        S.full_type=full_type;
        [ netmats ] = netmat_spectramt( data, S );
        sub_state_netmats.state{k}=netmats;
    end

    S=[];
    S.netmats=sub_state_netmats.global;
    S.type=type;
    S.full_type=full_type;
    [ netmats ] = netmat_spectramt( data, S );
    sub_state_netmats.global=netmats;
    
    state_netmats{ss}=sub_state_netmats;
end;