class PQueueNode:
    Next = None;
    data = None;
    def __init__(self, data):
        self.data = data
        self.Next = None

class PQueue:
    Head = None
    def __init__(self):
        self.Head = None
    def insert(self, data):
        if self.Head is None:
            self.Head = PQueueNode(data)
            return
        nNode = PQueueNode(data)
        if (data < self.Head.data):
            nNode.Next = self.Head
            self.Head = nNode
            return
        p = self.Head

        while(p.Next != None):
            if data > p.data and p.Next.data > data:
                nNode.Next = p.Next
                p.Next = nNode
                return
            p=p.Next
        p.Next = nNode
    def extract_min(self):
        res = self.Head
        self.Head = self.Head.Next
        res.Next = None
        return res
    def isempty(self):
        return (self.Head == None)