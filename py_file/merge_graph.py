import os
from PIL import Image

base_path = "/Users/hazel/Desktop/lab/first_task"
folder_2 = os.path.join(base_path, "sig_2sub_graphs")
folder_3 = os.path.join(base_path, "sig_3sub_graphs")
folder_4 = os.path.join(base_path, "sig_4sub_graphs")
output_folder = os.path.join(base_path, "sig_2-4subs_graphs")

os.makedirs(output_folder, exist_ok=True)

# Loop through files in sig_2sub_graphs
for filename in os.listdir(folder_2):
    if not filename.endswith(".png"):
        continue

    file_2 = os.path.join(folder_2, filename)
    file_3 = os.path.join(folder_3, filename)
    file_4 = os.path.join(folder_4, filename)

    # Check if the same file exists in all folders
    if os.path.exists(file_3) and os.path.exists(file_4):
        img2 = Image.open(file_2)
        img3 = Image.open(file_3)
        img4 = Image.open(file_4)
        
        width = max(img2.width, img3.width, img4.width)
        total_height = img2.height + img3.height + img4.height
        combined_img = Image.new('RGB', (width, total_height), (255, 255, 255))

        y_offset = 0
        for img in [img2, img3, img4]:
            combined_img.paste(img, (0, y_offset))
            y_offset += img.height

        combined_img.save(os.path.join(output_folder, filename))