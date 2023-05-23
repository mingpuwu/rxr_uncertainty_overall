
def GetReactionNum():
    fd = open('reactionnum')
    num = fd.read()
    return int(num)