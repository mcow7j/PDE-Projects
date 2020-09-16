function [u,v]=turdomlin(a,b,gamma,N,M,L,k1,uit,vit,Bxu,Byu,Bxv,Byv)

du=1;
dv=10;

    QoldU=QoldmakerUderiv(uit,vit,gamma,N,M);
    QU=QmakerUlin(uit,vit,QoldU,gamma,a,N,M);
    AU=Acalculator1d7(uit,L,Bxu,Byu);
    fu=(-uit-0.5*k1*du*AU-k1*QU);
    u = MultigridVL1timedbound8(uit,QoldU,fu,L,k1,1,Bxu,Byu);
    
   
    QoldV=QoldmakerVderiv(u,gamma,N,M);
    QV=QmakerVlin(u,vit,QoldV,gamma,b,N,M);
    AV=Acalculator1d7(vit,L,Bxv,Byv);
    fv=(-vit-0.5*k1*dv*AV-k1*QV);
    v =MultigridVL1timedbound8(vit,QoldV,fv,L,k1,10,Bxv,Byv);
    
    
   
end