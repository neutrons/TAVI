"""Tests for tavi.frontend.widget.column_combo."""

from tavi.frontend.widget.column_combo import set_combo_options, set_combo_text, typed_column_combo


def test_typed_column_combo_is_editable(qtbot):
    combo = typed_column_combo()
    qtbot.addWidget(combo)

    assert combo.isEditable()


def test_typed_column_combo_shows_its_default(qtbot):
    combo = typed_column_combo("detector")
    qtbot.addWidget(combo)

    assert combo.currentText() == "detector"


def test_typed_column_combo_keeps_typed_names_out_of_the_list(qtbot):
    combo = typed_column_combo()
    qtbot.addWidget(combo)
    combo.setEditText("made_up")

    assert combo.count() == 0


def test_set_combo_text_keeps_the_index_in_step(qtbot):
    combo = typed_column_combo()
    qtbot.addWidget(combo)
    combo.addItems(["s2", "detector"])

    set_combo_text(combo, "detector")

    assert combo.currentIndex() == 1
    assert combo.currentText() == "detector"


def test_set_combo_text_accepts_a_name_not_in_the_list(qtbot):
    combo = typed_column_combo()
    qtbot.addWidget(combo)
    combo.addItems(["s2"])

    set_combo_text(combo, "made_up")

    assert combo.currentIndex() == -1
    assert combo.currentText() == "made_up"


def test_set_combo_text_emits_nothing(qtbot):
    combo = typed_column_combo()
    qtbot.addWidget(combo)
    combo.addItems(["s2", "detector"])
    received = []
    combo.currentIndexChanged.connect(received.append)

    set_combo_text(combo, "detector")

    assert received == []


def test_set_combo_options_replaces_the_list(qtbot):
    combo = typed_column_combo()
    qtbot.addWidget(combo)
    combo.addItems(["s2"])

    set_combo_options(combo, ["en", "monitor"])

    assert [combo.itemText(i) for i in range(combo.count())] == ["en", "monitor"]


def test_set_combo_options_keeps_the_shown_name(qtbot):
    combo = typed_column_combo("detector")
    qtbot.addWidget(combo)

    set_combo_options(combo, ["s2", "detector"])

    assert combo.currentText() == "detector"


def test_set_combo_options_emits_nothing(qtbot):
    combo = typed_column_combo("detector")
    qtbot.addWidget(combo)
    received = []
    combo.currentIndexChanged.connect(received.append)

    set_combo_options(combo, ["s2", "detector"])

    assert received == []
