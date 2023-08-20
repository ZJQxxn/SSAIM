"
Description:  
    SPAV algorithm for smooth isotoic models.

Programmer: 
    Jiaqi Zhang

Reference:
    ``A smoothed monotonic regression via L2 regularization''
    -- https://link.springer.com/article/10.1007/s10115-018-1201-2

Functions:
    init: Initialization;
    gt: Compare operator '>';
    lt: Compare operator '<';
    blockValue: Calculate A and Y;
    tridiagonalSolver: Solve a linear equation with tridiagonal matrix;
    monotonicityTest: Test monotonicity for estimated responses;
    updateBlock: Update the block partitions;
    updateEstimateY: Updated estimated responses;

    **fitting**: Main function for SPAV algorithm.
"

#TODO: The correctness

init<-function(xVec,yVec,smoothPar=0.1){
    "
    Initialization
    "
    # Reordr the predictor vector.
    reorder=sort(xVec,index.return=TRUE)$ix
    estimateY=yVec
    sampNum=length(xVec)
    blocks=list()
    for (i in 1:sampNum){
        blocks[i]=c(i)
    }
    S=c() # A set.
    smoothVar=smoothPar
    monotone=FALSE
    return (list(reorder=reorder,estimateY=estimateY,sampNum=sampNum,blocks=blocks,S=S))
}

gt<-function(par1,par2){
    "
    Given two number 'par1' and 'par2', determine whether par1 > par2. 
    "
    if (par1 > par2){
        return (1)
    }
    else{
        return (0)
    }
}

lt<-function(par1,par2){
    "
    Given two number 'par1' and 'par2', determine whether par1 < par2. 
    "
    if (par1 < par2){
        return (1)
    }
    else{
        return (0)
    }
}

blockValue<-function(xVec,yVec,reorder,blocks,smoothPar){
    "
    Calculate matrix A and Y for each block.
    "
    sampNum=length(yVec)
    blkNum=length(blocks)
    #print(blkNum)
    A=matrix(0,blkNum,blkNum)
    Y=matrix(0,1,blkNum)
    for (index in 1:blkNum){
        # delta / block size
        coeff=smoothPar/length(blocks[index])
        # Y is the mean of each block
        #print(blocks[[index]])
        Y[index] = mean(yVec[reorder][blocks[[index]]])
        # A is a tridiagonal matrix
        A[index,index]=1+coeff*gt(blocks[[index]][1],1)+coeff*lt(blocks[[index]][length(blocks[[index]])],sampNum)
        if ((index+1)<=blkNum){
            A[index,index+1]=-coeff*lt(blocks[[index]][length(blocks[[index]])],sampNum)
        }
        if (index>1){
            A[index,index-1]=-coeff*gt(blocks[[index]][1],1)
        }
    }
    return (list(A=A,Y=Y))
}

tridiagonalSolver<-function(A,Y,blkNum){
    "
    Thomas algorithm to solve a linear system with a tridiagonal matrix.

    Reference:
        ``Wikipedia: Tridiagonal matrixc algorithm''
        -- https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    "        
    x=matrix(0,1,blkNum)
    estimate_coef=matrix(0,1,blkNum-1)
    estimate_response=matrix(0,1,blkNum)
    # Forward step
    if(blkNum>1){
        estimate_coef[1]=A[1,2]/A[1,1]
        estimate_response[1]=Y[1]/A[1,1]
    }else{
        estimate_coef[1]=A[1,1]/A[1,1]
        estimate_response[1]=Y[1]/A[1,1]
    }
    for (index in 2:blkNum){
        # Element number of estimate_response is more than that of estimate_coef
        if (index==blkNum){
            estimate_response[index] = (Y[index] - A[index, index - 1] * estimate_response[index - 1])/(A[index, index] - A[index, index - 1] * estimate_coef[index - 1])
        }
        else{
            estimate_coef[index]=A[index,index+1]/(A[index,index]-A[index,index-1]*estimate_coef[index-1])
            estimate_response[index]=(Y[index]-A[index,index-1]*estimate_response[index-1])/(A[index,index]-A[index,index-1]*estimate_coef[index-1])
        }
    }
    # Backward step
    x[length(x)]=estimate_response[length(estimate_response)]
    for (i in (blkNum - 1):1){
        x[i]=estimate_response[i]-estimate_coef[i]*x[i+1]
    }
    return (x)
}

monotonicityTest<-function(S,estimateY,reorder,monotone){
    "
    Test the monotonicity of estimated responses.
    "
    sampNum=length(estimateY)
    monotone=TRUE
    tempVec=estimateY[reorder]
    for (index in 1:(sampNum-1)){
        if (! index %in% S){
            if(tempVec[index] >= tempVec[index+1]){
                S=union(S,index)
            }
        }
        if (tempVec[index]>tempVec[index+1]){
            monotone=FALSE
        }
    }
    return (list(S=S,monotone=monotone))
}

updateBlock<-function(S,reorder,sampNum){
    "
    Update block partitions.
    "
    # Clear current block partitions.
    #print('S:')
    #print(S)
    blocks=list()
    tempBlk=c()
    index=1
    count=1
    while (index <= sampNum){
        p=reorder[index]
        tempBlk=c(tempBlk,p)
        #while (length(S)>0 && p %in% S){
        while (p %in% S){
            index = index+1
            if (index > sampNum){
                break
            }
            p = reorder[index]
            tempBlk=c(tempBlk,p)
        }
        blocks[[count]]=tempBlk
        count=count+1
        tempBlk=c()
        index=index+1
    }
    return (blocks)
}

updateEstimateY<-function(blkValue,blocks,estimateY){
    "
    Update the estimated response.
    "
    for (blkIndex in 1:length(blocks)){
        for (each in blocks[blkIndex]){
            estimateY[each]=blkValue[blkIndex]
        }
    }
    return (estimateY)
}

SPAV<-function(xVec,yVec,smoothPar,epoch=0){
    "
    Main function for SPAV algorithm.
    "
    # Initialization
    if (epoch==0){
        epoch=length(xVec)
    }
    initialization=init(xVec,yVec,smoothPar)
    reorder=initialization$reorder
    estimateY=initialization$estimateY
    sampNum=initialization$sampNum
    blocks=initialization$blocks
    S=initialization$S
    monotone=FALSE
    sampNum=length(yVec)
    rm(initialization)
    # Solve SSAIM
    curEpoch=1
    while (!monotone){
        blkNum=length(blocks)
        if(blkNum==1){
            break
        }
        blk= blockValue(xVec,yVec,reorder,blocks,smoothPar)
        A=blk$A
        Y=blk$Y
        rm(blk)
        blkValue = tridiagonalSolver(A,Y,blkNum)  
        estimateY=updateEstimateY(blkValue,blocks,estimateY)
        test=monotonicityTest(S,estimateY,reorder,monotone)
        S=test$S
        monotone=test$monotone
        rm(test)
        if (monotone){
            break
        }
        blocks=updateBlock(S,reorder,sampNum)

        curEpoch=curEpoch+1
        if (curEpoch>=epoch){
            #print('Reach at the predefined iteration number.')
            break
        }
    }
    #print(sprintf('Finish SPAV in %d steps!',curEpoch))
    return (estimateY)
}
