# PyTorch network

import torch.nn as nn
import torch.nn.functional as F

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = nn.Conv2d(in_channels = 1, out_channels = 8, kernel_size = 16, stride = (2,2), padding = 7)
        
        self.conv2 = nn.Conv2d(in_channels = 8, out_channels = 8, kernel_size = 8, padding = 1)
        
        self.pool1 = nn.MaxPool2d(kernel_size = (2,2), stride = (2,2))
        
        self.conv3 = nn.Conv2d(in_channels = 8, out_channels = 8, kernel_size = 7, padding = 3) # input 20x20x8
        
        self.pool2 = nn.MaxPool2d(kernel_size = (2,2), stride = (2,2))
        
        self.conv4 = nn.Conv2d(in_channels = 8, out_channels = 8, kernel_size = 4, stride = (2,2), padding = 1)
        
        self.avepool = nn.AvgPool2d(kernel_size = (2,2), stride = (1,1))
        
        self.fc1 = nn.Linear(8*4*4, 200)
        
        self.fc2 = nn.Linear(200, 200)
        
        self.fin = nn.Linear(200, 1)
        
        self.sig = nn.Sigmoid()

    def forward(self, x):
        x = F.relu(self.conv1(x))

        x = self.pool1(F.relu(self.conv2(x)))

        x = self.pool2(F.relu(self.conv3(x)))

        x = self.avepool(F.relu(self.conv4(x)))
        #print(x.size())

        x = x.view(x.size(0), 8*4*4)

        x = F.relu(self.fc1(x))

        x = F.relu(self.fc2(x))

        x = self.fin(x)
        x = self.sig(x)

        return x

