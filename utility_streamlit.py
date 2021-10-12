import os, json
# @st.cache(allow_output_mutation=True)
def load_global_recovered_values():
    print('[load_global_recovered_values] for cache')
    # global_recovered_values = dict()
    with open(f'{os.path.dirname(__file__)}/global_recovered_values.json', 'r') as f:
        global_recovered_values = json.load(f)
    return global_recovered_values
def make_recording_widget(f, global_recovered_values):
    # https://discuss.streamlit.io/t/how-to-save-all-widget-state/1764/2
    """Return a function that wraps a streamlit widget and records the widget's values to a global dictionary.
    """
    def wrapper(label, *args, **kwargs):
        # print('The label:', label)
        widget_value = f(label, *args, **kwargs)
        global_recovered_values[label] = widget_value
        return widget_value
    return wrapper
