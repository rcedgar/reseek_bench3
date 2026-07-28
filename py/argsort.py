def argsort(arr, reverse=False):
    return sorted(range(len(arr)), key=lambda i: arr[i], reverse=reverse)
