function [met]=QuadMutNetExMetricsStruct(sol,mx,subNetworkSolution,c,k,np)

subNetworkSolution = np * subNetworkSolution;

solMx=full(mx(:,sol));
    [sCnt,gCnt]=size(solMx);
    covMx=sum(solMx,2);
    mCnt=sum(covMx);
    covTot=length(find(covMx>0));
    covR=covTot/sCnt;
    covOvr=mCnt-covTot;
    excessR=length(find(covMx>1))./covTot;
    lObj=covTot-covOvr;
    
    G=solMx;
    Q=transpose(G) * G;
    Q=Q-subNetworkSolution; % NET
	H=.5*(Q+transpose(Q));
	H=.5*(k+1)*H;
	f=(-((k+3)/2) *ones(1,sCnt))*G;
    f=f+c.*ones(size(f));
    x=ones(1,gCnt);
    qObj=sCnt+x*H*(x')+f*(x');
    met.qObj=qObj;
    met.lObj=lObj;
    met.covTot=covTot;
    met.covOvr=covOvr;
    met.covR=covR;
    met.excessR=excessR;
    met.mCnt=mCnt;
    met.sCnt=sCnt;
    met.gCnt=gCnt;
	    
