
#reference code from ds9_reader
from unittest import skip


#source_ids = []
#for epoch in new_epochs:
#    source_ids.append(new_epochs[epoch][0][:9])


image_list = []

for sid in source_ids:
    if sid != files:
        skip
    if sid == files:
        image_list.append(sid)
        if len(image_list) == 2:
            if sid[10:12] == 'w2':
                image_list.append(image_list.pop(sid))
            if sid[10:12] == 'w1':
                image_list.pop(sid)
                image_list.insert(0, sid)

