const { create, all } = require('mathjs');
const {  performance } = require('perf_hooks');
const config = { }
const math = create(all, config)

function isPosDef(A) {
    return  math.eigs(A).values.every(value => value > 0)
}

function changeSign(A) {
    return A.map(function (value) {
        return value.map(function (value) {
            return   value * (-1)
        })
    })
}
let arraysMatch = function (arr1, arr2) {

    // Check if the arrays are the same length
    if (arr1.length !== arr2.length) return false;
    // Check if all items exist and are in the same order
    for (let i = 0; i < arr1.length; i++) {
        for(let j = 0; j< arr1.length; j++) {
            if (arr1[i][j] !== arr2[i][j]) return false;
        }
    }
    return true;

};


function gradientDescent (A, b, x,e) {
    let transposeA = math.transpose(A)
    console.log(arraysMatch(A,transposeA))
    if(!isPosDef(A) || !arraysMatch(A, transposeA)) {
        throw Error("Matrix A needs to be symmetric positive definite")
    } else {
        //calculates inverse of matrix
        let invA = math.inv(A);
        //calculates the exact min
        // let signInvA = changeSign(invA);
        // let invAb = math.multiply(signInvA,b);
        // console.log('the real min: ',invAb);
        //let Ax = math.multiply(A,x);
        let signA = changeSign(A);
        let signAx = math.multiply(signA,x);
        let minusB = b.map(value => value * (-1));
        // h(xk) = -f'(xk);
        let gradient = math.subtract(minusB,signAx);
        console.log('h(xk): ',gradient)
        let k = 0;
        //while norm of derivative is bigger than epsilon, build xk+1;
        //minimization of function with recurrent relation
        while(math.norm(gradient) > e) {
            //p is hk
            let p = gradient;
            console.log('hk: ',p)
            // the q is A*hk
            let q = math.multiply(signA,p);
            console.log('A*hk: ', q)
            //alphaK = (f'(xk),h(k)/(Ahk,hk);
            let alpha = math.dotDivide(math.dot(p,gradient),math.dot(q,p));
            console.log('alpha ', alpha)
            let z = math.dotMultiply(alpha,p);
            //minimization of function with  recurrent relation
            x = math.subtract(x,z);
            console.log('x: ', x)
            let h = math.dotMultiply(alpha,q);
            //To avoid multiplying by A twice per iteration,
            //it x = x = alpha*A*hk implies hk = hk - alpha * A * hk
            //it can be gradient = -b - Axk
            gradient = math.subtract(gradient,h);
            console.log('gradient ',gradient)
            k += 1;
            console.log('iterations number is: ',k)
        }
        console.log('x ::',x)
        return x
    }
}

let t0 = performance.now()
gradientDescent([[3, 2], [2, 3]],[-2,7],[0,0],1e-10)
let t1 = performance.now()
console.log("Call to gradientDescent took " + (t1 - t0) + " milliseconds.")

