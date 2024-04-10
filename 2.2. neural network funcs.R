#######################################################
#Functions to fit a neural network
#######################################################

#the neuralnet package is not friendly with user-defined loss function
#the keras package does not run well on the server
#we need to write our own feed forward neural network functions

#------------------------
#Layer functions
#------------------------

#' Activation functions
#' 
#' @description Default activation functions
#' @param x a vector of input data
#' @param type type of activation function. Options are "identity", "relu", "tanh", "sigmoid", "softmax"
#' @param derivative whether the derivative is returned. Default is FALSE, meaning that the value of the function is returned.
#' @return value of the activation function or its derivative at data x
t_activate <- function(x, type = c("identity", "relu", "tanh",
                                   "sigmoid", "softmax"), 
                       derivative = FALSE) {
  type <- match.arg(type)
  
  if (type == "identity") {
    if (derivative) {
      ret <- rep(1, length(x))
    } else {
      ret <- x
    }
  } else if (type == "relu") {
    if (derivative) {
      ret <- as.numeric(x > 0)
    } else {
      ret <- x
      ret[x <= 0] <- 0
    }
  } else if (type == "tanh") {
    ret <- tanh(x) #(exp(x) - exp(-x))/(exp(x) + exp(-x))
    if (derivative) {
      ret <- 1 - ret^2
    } 
  } else if (type == "sigmoid") {
    ret <- 1/(1+exp(-x))
    if (derivative) {
      ret <- ret*(1-ret)
    }
  } else if (type == "softmax") {
    tmp <- exp(x)
    ret <- tmp/sum(tmp)
    if (derivative) {
      ret <- sweep((diag(length(x)) - ret), MARGIN = 2, ret, `*`) 
    }
  }
  ret
}

#' Loss functions
#' 
#' @description Default loss functions
#' @param y true values of response variable
#' @param ypred the predicted value of response variable
#' @param type type of loss function. Options are "mse", "ce" (cross entropy) and "bernoulli" (log likelihood of independent bernoulli distribution)
#' @param derivative whether the derivative is returned. Default is FALSE, meaning that the value of the function is returned.
#' @param clip ypred may have irregular values that is not defined for the loss. For example, "ce" and "bernoulli" requires ypred to be between 0 and 1. Clip constrains ypred to be between (clip, 1-clip).
#' @return value of the loss function or its derivative
t_loss <- function(y, ypred, 
                   type = c("mse", "ce", "bernoulli"), 
                   derivative = FALSE, clip = 1e-08) {
  
  type <- match.arg(type) #type of loss function
  
  if (type == "mse") {
    if (derivative) {
      ret <- -2*(y - ypred)
    } else {
      ret <- (y - ypred)^2 #minus log likelihood of normal distribution
    }
  } else if (type == "ce") {
    #clipping to ensure numerical stability
    ypred <- sapply(1:length(ypred), function(i) {max(clip, ypred[i])})
    ypred <- sapply(1:length(ypred), function(i) {min(ypred[i], 1-clip)})
    if (derivative) {
      ret <- -y/ypred
    } else {
      ret <- -sum(y*log(ypred)) #minus log likelihood of multinomial/categorical distribution
    }
  } else if (type == "bernoulli") {
    #clipping to ensure numerical stability
    ypred <- sapply(1:length(ypred), function(i) {max(clip, ypred[i])})
    ypred <- sapply(1:length(ypred), function(i) {min(ypred[i], 1-clip)})
    if (derivative) {
      ret <- -(y/ypred - (1-y)/(1-ypred))
    } else {
      ret <- -sum(y*log(ypred) + (1-y)*log(1-ypred)) #minus log likelihood of indep bernoulli dists
    }
  }
  
  ret #return
}

#' Add intercept to a input vector
#' 
#' @description This function adds intercept (bias) to the output of a layer to input to the next
#' @param vec a vector
#' @return an augmented vector with intercept (1) added on top
t_addBias <- function(vec) {
  c(1, vec)
}

#------------------------
#Basic training functions
#------------------------

#' Move input forward one layer
#' 
#' @description This function moves forward the input from one layer to the next
#' @param input a vector of one observation input (x)
#' @param actFUN activation function of the current layer
#' @param weight weight matrix connecting current layer to the next
#' @return the output from the current layer, which can be used as the input to the next layer
t_forward <- function(input, actFUN, weight) {
  
  #input, actFUN of layer i
  #weight connecting layer i and layer i+1
  
  output <- actFUN(input)
  if (dim(weight)[2] > length(input)) {
    output <- weight %*% as.matrix(t_addBias(output))
  } else {
    output <- weight %*% as.matrix(output)
  }
  
  as.vector(output) #return
}

#' Backpropagation in one layer 
#' 
#' @description This function conducts a backpropagation of loss derivative from one layer to the next
#' @param delta the loss derivative according to the current layer
#' @param weight the weight that acconect the current layer to the previous layer
#' @param input the input at the previous layer
#' @param actFUN the activation function at the previous layer
#' @param actFUN.derivative the activation function at the previous layer
#' @return an object that contains the loss derivative according to the previous layer, and the gradient for the weight connecting current layer to the previous layer
t_backward <- function(delta, weight, input, 
                       actFUN, actFUN.derivative) {
  
  #delta is the derivative of loss wrt input layer i+1
  #weight i connecting layer i and i+1
  #input, actFUN and actFUN.derivative of layer i
  
  activated <- actFUN(input)
  if (dim(weight)[2] > length(activated)) {
    activated <- t_addBias(activated)
  }
  
  derivative <- actFUN.derivative(input)
  if (dim(weight)[2] > length(derivative)) {
    derivative <- c(0, derivative)
  }
  
  gradient <- as.matrix(delta) %*% t(as.matrix(activated)) 
  tmp <- sweep(weight, MARGIN = 2, derivative, FUN = `*`)
  delta <- sweep(tmp, MARGIN = 1, delta, FUN = `*`)
  delta <- colSums(delta) #delta now contains bias derivative
  
  list(delta = delta[2:length(delta)], gradient = gradient)
}

#------------------------
#Feed forward neural network object
#------------------------

#' Create a feed forward neural network
#' 
#' @description This function creates an object that contains the structure of a feed forward neural network 
#' @param layers a numeric vector that contains the number of neurons per input layer + hidden layer + output layer
#' @param actFUN a list of activation functions for each layer. If a character vector is provided, then the corresponding types of default activation functions provided in t_activate() is used
#' @param actFUN.derivative a list of corresponding derivative of activation functions. If a character vector is provided, then the corresponding types of derivatives in the t_activate() is used 
#' @param lossFUN the loss function. If a character string is provided, then the corresponding type of loss function provided in default function t_loss will be used
#' @param lossFUN.derivative the derivative of the loss function. If a character string is provided, then the corresponding type of derivatives provided in default function t_loss will be used
#' @param learning_rates the learning rate in each layer
#' @param weights the weights connecting layers in the network. The default is NULL. That is, the weights are initialized randomly.
#' @return an object of class "ffnn" that contains relevant information on the structure of the neural network.
ffnn <- function(layers, actFUN, actFUN.derivative = actFUN, 
                 lossFUN, lossFUN.derivative = lossFUN, 
                 learning_rates = NULL, weights = NULL) {
  
  #weights
  if (is.null(weights)) {
    weights <- list()
    for (i in 1:(length(layers)-1)) {
      #default has bias
      weights[[i]] <- (0.2*matrix(rnorm((layers[i]+1)*layers[i+1]), 
                                  nrow = layers[i+1])-0.1)/10000
    } 
  }
  
  #activation function
  if (is.vector(actFUN.derivative)) {
    actFUN.derivatives <- list()
    for (i in 1:length(layers)) {
      actFUN.derivatives[[i]] <- t_activate
      formals(actFUN.derivatives[[i]])[c("type", "derivative")] <- list(type = actFUN.derivative[i], 
                                                                        derivative = TRUE)
    }
    actFUN.derivative <- actFUN.derivatives
  }
  
  if (is.vector(actFUN)) {
    actFUNs <- list()
    for (i in 1:length(layers)) {
      actFUNs[[i]] <- t_activate
      formals(actFUNs[[i]])$type <- actFUN[i]
    }
    actFUN <- actFUNs
  }
  
  #loss function
  if (is.character(lossFUN.derivative)) {
    lossFUN.derivatives <- t_loss
    formals(lossFUN.derivatives)$type <- lossFUN.derivative
    formals(lossFUN.derivatives)$derivative <- TRUE
    lossFUN.derivative <- lossFUN.derivatives
  }
  
  if (is.character(lossFUN)) {
    lossFUNs <- t_loss
    formals(lossFUNs)$type <- lossFUN
    lossFUN <- lossFUNs
  } 
  
  #learning rates
  if (is.null(learning_rates)) {
    learning_rates <- rep(0.01, length(weights))
  }

  #return
  ret <- list(layers = layers, weights = weights, 
              actFUN = actFUN, actFUN.derivative = actFUN.derivative, 
              lossFUN = lossFUN, lossFUN.derivative = lossFUN.derivative, 
              learning_rates = learning_rates)
  class(ret) <- "ffnn"
  ret
}

#------------------------
#Training functions
#------------------------

#' Gradients of the neural network weights
#' 
#' @description This function calculates gradients for all weights of the feed forward neural network based on one input observation
#' @param obj an object of class "ffnn"
#' @param input a vector of one input observation
#' @param y the corresponding response value for the input observation
#' @return an object that contains the predicted value, gradients (the set of gradients for each weight matrix), and the loss  
gradient.ffnn <- function(obj, input, y) {
  
  #feedforward through all layers of the neural network
  layers <- list()
  layers[[1]] <- input
  for (i in 1:(length(obj$layers)-1)) {
    layers[[i+1]] <- t_forward(layers[[i]], 
                               actFUN = obj$actFUN[[i]], 
                               weight = obj$weights[[i]])
  }
  
  #get the prediction
  ypred <- obj$actFUN[[length(obj$layers)]](layers[[length(obj$layers)]])
  ypred <- as.vector(ypred)
  
  #delta from loss function
  delta_loss <- obj$lossFUN.derivative(y = y, ypred = ypred)
  loss <- obj$lossFUN(y = y, ypred = ypred)
  
  #delta from last layer
  delta <- obj$actFUN.derivative[[length(obj$layers)]](layers[[length(obj$layers)]])
  if (is.matrix(delta)) {
    delta <- as.vector(delta %*% delta_loss)
  } else {
    delta <- delta * delta_loss
  }
  
  #backpropagation
  gradients <- list()
  for (i in length(obj$weights):1) {
    tmp <- t_backward(delta = delta, weight = obj$weights[[i]], 
                      input = layers[[i]], actFUN = obj$actFUN[[i]], 
                      actFUN.derivative = obj$actFUN.derivative[[i]])
    delta <- tmp$delta
    gradients[[i]] <- tmp$gradient
  }
  
  list(ypred = ypred, gradients = gradients, loss = loss)
}

#' Gradient accumulation 
#' 
#' @description This function sums two sets of gradients 
#' @param gradients1 the existing set of gradients (or the set of gradients accumulated so far) 
#' @param gradients2 the new set of gradients to be added to gradients1
#' @param factor the factor by which gradients2 is multiplied before being added to gradients1
#' @return a new set of gradients that equals gradients1 + factor*gradients2
accumulate_weight.ffnn <- function(gradients1, gradients2, factor = 1) {
  
  for (i in 1:length(gradients1)) {
    gradients1[[i]] <- gradients1[[i]] + factor*gradients2[[i]] 
  }
  
  gradients1
}

#' Weight updating
#' 
#' @description This function updates the weights of the neural network using sets of gradients and learning rates
#' @param obj an "ffnn" (feed forward neural network) object
#' @param gradients the set of gradients to be used to update the weights of the ffnn obj
#' @param learning_rates the set of learning rates for each weight in obj. The default is NULL, which means the learning rates stored in the ffnn obj will be used.
#' @return the ffnn obj with weights updated.
update.ffnn <- function(obj, gradients, 
                        learning_rates = NULL) {
  
  #possibility for a prespecified learning rates, other than object's one
  if (is.null(learning_rates)) {
    learning_rates <- obj$learning_rates 
  }
  
  for (i in 1:length(obj$weights)) {
    obj$weights[[i]] <- obj$weights[[i]] - learning_rates[i]*gradients[[i]] 
  }
  
  obj
}

#' Training the feed forward neural network
#' 
#' @description This function trains a feed forward neural network using a training data set
#' @param obj the ffnn object
#' @param input the training input matrix (x) of size nxp, where n is the number of observations, p is number of independent variables
#' @param y the training response values (y)
#' @param batchsize the training batchsize. Default is 25
#' @param epoch the number of training epochs (iteration). Default is 100
#' @param verbose whether the training process is printed. 
#' @param decay the rates by which the learning rates decrease in each epoch (iteration)
#' @return the ffnn object with weights trained on the training data.
train.ffnn <- function(obj, input, y,  
                       batchsize = 25, epoch = 100, 
                       verbose = TRUE, decay = 0) {
  
  n.input <- dim(input)[1]
  n.batch <- ceiling(n.input/batchsize)
  y <- as.matrix(y)
  
  for (i in 1:epoch) {
    
    #randomize order of data for each epoch
    samp <- sample(1:n.input, n.input, replace = FALSE)
    y.train <- as.matrix(y[samp,])
    x.train <- input[samp,]
    
    #calculate learning rates for this epoch
    learning_rates <- obj$learning_rates/(1+decay*i)
    
    #initialize loss for this epoch
    loss <- c()

    for (j in 1:n.batch) {
      
      #initialize gradients
      gradients <- as.list(rep(0, length(obj$weights)))
      
      #calculate gradient for each input in the batch
      for (k in ((batchsize*(j-1)+1):(min(batchsize*j, n.input)))) {
        
        #current input
        input.tmp <- x.train[k,]
        y.tmp <- y.train[k,]
        
        #calculate gradients using backpropagation
        tmp <- gradient.ffnn(obj, input = input.tmp, y = y.tmp)
        
        #update the gradients
        for (l in 1:length(gradients)) {
          gradients[[l]] <- gradients[[l]] + tmp$gradients[[l]]
          loss <- c(loss, tmp$loss)
        }
      }
      
      #update weights per batch
      obj <- update.ffnn(obj, gradients = gradients, learning_rates = learning_rates)
    }
    
    #print progress
    if (verbose) {
      cat("Epoch ", i, "/", epoch, "loss = ", mean(loss), "\n")
    }
  }

  obj #return obj with updated weights
}

#' Prediction using feed forward neural network
#' 
#' @description This function uses a feed forward neural network to predict responses based on input
#' @param obj the ffnn object
#' @param input the input (x) matrix of size nxp, where n is the number of observations, p is number of independent variables
#' @return the predicted responses
predict.ffnn <- function(obj, input) {
  
  n.input <- dim(input)[1]
  output <- c()
  
  for (i in 1:n.input) {
    input_tmp <- input[i,]
    for (i in 1:(length(obj$weights))) {
      input_tmp <- t_forward(input = input_tmp, 
                             actFUN = obj$actFUN[[i]], 
                             weight = obj$weights[[i]])
    }
    output_tmp <- obj$actFUN[[length(obj$layers)]](input_tmp)
    output <- rbind(output, output_tmp)
  }
  
  rownames(output) <- 1:n.input
  output #return
}

#------------------------
#Ensemble of networks
#------------------------

#' Initialize an ensemble of neural networks
#' 
#' @description This function initializes an ensemble of feed forward neural networks of the same architecture 
#' @param ffnn_args parameters to be inputted in the ffnn() function. That is, the mutual structure of the neural networks 
#' @param rep the number of neural networks in the ensemble
#' @return an object of class "ensemble" that contains the list of ffnn objects with similar structure.
ensemble_ffnn <- function(ffnn_args, rep = 10) {
  
  ret <- list()
  for (i in 1:rep) {
    ret[[i]] <- do.call(ffnn, ffnn_args)
  }
  
  class(ret) <- "ensemble"
  ret
}

#' Train an ensemble of neural networks
#' 
#' @description This function trains an ensemble of neural networks
#' @param obj an object of class "ensemble", i.e., the output of the ensemble_ffnn() function 
#' @param input the training input matrix (x) of size nxp, where n is the number of observations, p is number of independent variables
#' @param y the training response values (y)
#' @param verbose whether the training process (for each neural network) is printed
#' @param ... other parameters to be passed to the train.ffnn() function
#' @return an ensemble of networks with weights trained on training data
train.ensemble <- function(obj, input, y, verbose = FALSE, ...) {
  
  nffnn <- length(obj)
  for (i in 1:nffnn) {
    obj[[i]] <- train.ffnn(obj[[i]], input = input, y = y, verbose = verbose, ...)
  }
  
  obj
}

#' Prediction from an ensemble of networks
#' 
#' @description This function uses an ensemble of neural networks to make prediction
#' @param obj an object of class "ensemble" that contains the ensemble of networks
#' @param input the input (x) matrix of size nxp, where n is the number of observations, p is number of independent variables
#' @param type the type of prediction method. It can be either UCB (default), Thompson sampling, or the posterior mean.
#' @return the predicted responses
predict.ensemble <- function(obj, input, 
                             type = c("ucb", "thompson", "postmean")) {
  
  #let each model in the ensemble make a prediction
  output.list <- list()
  nffnn <- length(obj)
  for (i in 1:nffnn) {
    output.list[[i]] <- predict.ffnn(obj[[i]], input)
  }
  
  #output according to type of the prediction
  type <- match.arg(type)
  if (type == "thompson") {
    #thompson sampling
    model.idx <- sample(1:nffnn, 1) 
    output <- output.list[[model.idx]]
  } else {
    
    vec <- unlist(output.list, use.names = FALSE)
    DIM <- dim(output.list[[1]])
    n <- length(output.list)
    
    list.mean <- tapply(vec, rep(1:prod(DIM),times = n), mean)
    attr(list.mean, "dim") <- DIM
    
    if (type == "postmean") {
      #posterior mean
      output <- list.mean
    } else if (type == "ucb") {
      #upper confidence bound
      list.sd <- tapply(vec, rep(1:prod(DIM),times = n), sd)
      attr(list.sd, "dim") <- DIM
      output <- list.mean - list.sd #want minimum
    }
  }
  
  output
}

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(t_activate, t_addBias, t_backward, t_forward, t_loss,
     ffnn, ensemble_ffnn, train.ensemble, train.ffnn, 
     gradient.ffnn, accumulate_weight.ffnn, update.ffnn, 
     predict.ffnn, predict.ensemble, 
     file = "funcs/neuralnet.RData")

#######################################################
#Test if things work as expected
#######################################################

load("funcs/neuralnet.RData")

#---------------------------
#softmax classification
#---------------------------

set.seed(123456)
data <- iris[sample(nrow(iris)),] #iris data
y <- data[, "Species"] #response 
x <- data[, 1:4] #input

#scale to [0,1]
x <- as.matrix(apply(x, 2, function(x) (x-min(x))/(max(x) - min(x))))

#one hot encode classes of labels
y <- as.integer(y)
y <- cbind(ifelse(y == 1, 1, 0), 
           ifelse(y == 2, 1, 0), 
           ifelse(y == 3, 1, 0))

model <- ffnn(layers = c(4, 3, 3), actFUN = c("relu", "relu", "softmax"), 
              lossFUN = "ce", learning_rates = c(0.1, 0.1))
tmp <- train.ffnn(obj = model, input = x, y = y, epoch = 100, decay = 0.2)

#check some prediction
idx <- 3
predict.ffnn(tmp, t(as.matrix(x[idx,]))) #true label receives highest score 
y[idx,]

#---------------------------
#regression
#---------------------------

set.seed(123456)
data <- iris[sample(nrow(iris)),]
y <- data[, 4]
x <- data[, 1:3]

# scale to [0,1]
x <- as.matrix(apply(x, 2, function(x) (x-min(x))/(max(x) - min(x))))

model <- ffnn(layers = c(3, 3, 1), actFUN = c("relu", "relu", "identity"), 
              lossFUN = "mse", learning_rates = c(0.01, 0.01))
tmp <- train.ffnn(obj = model, input = x, y = y, epoch = 100, decay = 0.1)

#check some prediction
idx <- 1:10
predict.ffnn(tmp, x[idx,]) #ok prediction results
y[idx]

#compare with neural net
library(neuralnet)
tmp2 <- neuralnet(Petal.Width ~ Sepal.Width + Sepal.Length + Petal.Length, data = data)
predict(tmp2, newdata = x[idx,]) #bad results, probably need hyperparameter tuning
y[idx] 

#compare with linear regression
tmp3 <- lm(Petal.Width ~ Sepal.Width + Sepal.Length + Petal.Length, data = data)
predict(tmp3, newdata = data[idx,]) #good results here as well
y[idx]


