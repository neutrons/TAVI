# AGENTS.md

Guidance for AI coding agents working in this repo. TAVI (Triple-Axis data
VIsualization toolkit) is a PyQt/PySide desktop app for visualizing and
reducing triple-axis spectrometer data.

## Environment & commands

This repo uses **Pixi**, not bare `pip`/`python`. Always run through it:

```bash
pixi install                 # sync/create the environment
pixi shell                   # enter it interactively
pixi run tavi                # launch the application
pixi run test                # run the test suite (pytest, integration tests excluded)
pixi run pytest <path>       # run a subset directly (pytest task also works with args)
pixi run ruff check <path>   # lint
pixi run ruff format <path>  # format
pixi run build-docs          # build Sphinx docs to docs/_build/html
pixi run test-docs           # run doctest build
```

`pixi task list` shows everything else (`clean`, `audit-deps`, `conda-build`, ...).

Integration tests are excluded by default (`addopts = "-m 'not (integration)'"`
in `pyproject.toml`); run them explicitly with `pixi run pytest -m integration`.

## Code layout

```
src/tavi/
  frontend/
    presenter/   # Mediates between EventBroker and a view. Holds no domain data.
    view/        # PyQt widgets. Pure UI — no business logic, no direct model access.
    widget/      # Reusable widget pieces.
  backend/
    model/       # Domain state + event handlers (TaviProjectModel, PlotModel, ...).
  library/
    data/        # Pydantic domain models (Plot, Scan, RawScan, TaviData, UUID, ...).
    storage/     # Loaders and file store abstractions.
    tas/, resolution/, geometry/, component/, experiment/, fit/, plot/, ubalgorithm/
                 # Triple-axis science/math library, independent of the GUI.
  meta/
    event/       # EventBroker + typed Event classes (the app's nervous system).
    exception/, multithreading/, decorators/, logging
tests/
  unit/          # Mirrors src/ layout.
  integration/   # Marked `integration`, excluded from default test run.
docs/source/     # Sphinx docs — guides/ (how things work) and design/ (why).
```

## Architecture: event-driven MVP

This is **not** a direct-call MVP — presenters and models never hold live
references to each other's storage. Everything flows through the singleton
`EventBroker` (`tavi.meta.event.event_broker.EventBroker`):

- Every event is a pydantic `Event` subclass in `tavi/meta/event/type/`.
- `EventBroker().register(EventType, handler)` subscribes; `.publish(event)`
  dispatches synchronously to every subscriber, each getting a **deep copy**
  of the event (`model_copy(deep=True)`) — so an event can safely carry
  model-owned domain objects without callers being able to mutate shared state.
- Dispatch is synchronous and recursion-guarded (`call_depth_max`, currently
  5) — a handler publishing another event nests one level deeper. Read
  `docs/source/guides/event_broker.rst` before adding a new link in an
  existing chain; it's easy to blow the depth budget.
- **Presenters hold uuids, not domain objects.** A presenter may cache a
  `list[UUID]` between events, but never a live `Plot`/`Scan`/`RawScan`. Data
  a presenter needs to render travels *inside* the event that triggers it.
  Models (`TaviProjectModel`, `PlotModel`) are the single source of truth for
  domain data (`TaviData.raw_scans`, `TaviData.plots`, `PlotModel._last_plots`).
- Read `docs/source/design/frontend/visualization_flow.rst` for a worked,
  diagrammed example of a full event chain (tree selection → rendered plot),
  including the lighter-weight "switch the active plot" side-channel.

When adding a new interaction, prefer adding a new narrowly-scoped `Event`
type over overloading an existing one — see "Prefer domain-specific events"
in `docs/source/guides/event_broker.rst`.

## Testing conventions

- Tests mirror `src/` layout under `tests/unit/`.
- Pytest + `pytest-qt` (`qtbot` fixture) for widget tests.
- Construct domain objects (`Plot`, `RawScan`, `UUID`, ...) via small
  `make_*()` helper functions at the top of each test file rather than
  hand-rolling pydantic kwargs inline — follow the existing pattern in
  neighboring test files.
- When testing event-driven behavior, register a throwaway
  `received = []; EventBroker().register(SomeEvent, received.append)` and
  assert on what got published, rather than mocking the broker.
- `EventBroker` is a singleton, but `tests/conftest.py` has an autouse
  fixture (`_reset_Singletons`, via `neutrons_standard`) that resets all
  singletons — including the broker's subscriber registry — before every
  test, class, and module. You don't need to reset it yourself.

## Style notes specific to this repo

- Every public class/function gets a one-line (occasionally two-line)
  docstring stating what it does; ruff's `D` (pydocstyle) rules are enforced.
- Domain models are pydantic `BaseModel`s (`ConfigDict(arbitrary_types_allowed=True)`
  where needed for a custom `UUID` type) — prefer `model_copy(update={...})`
  over mutation when producing a modified copy.
- Line length 120 (ruff `line-length = 120`).
- Comments explain *why*, not *what* — this matches the file's own style
  throughout `src/`; don't add narration comments to straightforward code.
- Pre-commit runs `ruff check --fix`, `ruff-format`, and `codespell` on
  `src/` (tests/scripts/prototype are excluded from ruff via
  `.pre-commit-config.yaml`, but keep them clean anyway).

## Docs

Sphinx docs live in `docs/source/`:
- `guides/` — practical how-it-works references (event broker, loaders, etc.)
- `design/` — architectural rationale, sequence diagrams (Mermaid via
  `.. mermaid::` directive), key-decisions sections

When a change alters an event flow, a model's ownership boundary, or a
presenter/view contract that a `design/` doc describes, update that doc in
the same change — these docs go stale otherwise (this file's own accuracy
depends on the same discipline). Diagrams and prose in `design/*.rst` are the
canonical source of "why it's built this way"; don't let code and docs drift.
