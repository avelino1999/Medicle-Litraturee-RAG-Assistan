import google.generativeai as genai

genai.configure(api_key="AIzaSyA7OOTvvQKFWR1luszOeBJla1S7ZIiP-Wk")

for model in genai.list_models():
    print(model.name)
