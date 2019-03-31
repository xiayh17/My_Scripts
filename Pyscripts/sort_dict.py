def sort_dict_by_value(mydict, reverse=True):
	'''
	sort dict by value
	'''
	from collections import OrderedDict
	if isinstance(mydict, dict):
		order_dict = OrderedDict(sorted(mydict.items(), key=lambda x:x[1],reverse = reverse))
		return order_dict
	else:
		raise TypeError("input must be a dict")
