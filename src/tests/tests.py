from contextvars import copy_context
from dash._callback_context import context_value
from dash._utils import AttributeDict

# Import the names of callback functions you want to test
from src.pages.download import make_dl_table, generate_download


def test_make_dl_table():
    output = make_dl_table(url=1)
    assert output == "button 1: 1 & button 2: 0"


def test_generate_download():
    def run_callback():
        context_value.set(
            AttributeDict(
                **{"triggered_inputs": [{"prop_id": "btn-1-ctx-example.n_clicks"}]}
            )
        )
        return generate_download(1, 0, 0)

    ctx = copy_context()
    output = ctx.run(run_callback)
    assert output == f"You last clicked button with ID btn-1-ctx-example"
