

from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.event_interface import Event

import pytest
import unittest


class DummyEvent(Event):
    pass

class SmartyEvent(Event):
    pass


class TestEventBroker(unittest.TestCase):
    def setUp(self):
        self.event_broker : EventBroker = EventBroker()
        pass

    def tearDown(self):
        pass
    
    def test_register(self):
        def dummy_method(_):
            pass
        self.event_broker.register(DummyEvent, dummy_method)
        
        assert len(self.event_broker.registry) == 1
        assert len(self.event_broker.registry[DummyEvent]) == 1
        assert self.event_broker.registry[DummyEvent][0] == dummy_method
        
    def test_publish(self):
        assert len(self.event_broker.registry) == 0, "singleton not cleared before start of test"
        e = DummyEvent()
        method_called = False
         
        def dummy_method(event:DummyEvent):
            nonlocal method_called, e
            method_called = True
            assert type(event) is DummyEvent, "wrong event type forwarded to subscriber"
            assert e is not event, "deep copy of event failed"
            
            
        self.event_broker.register(DummyEvent, dummy_method)
        self.event_broker.publish(e)
        
        assert method_called, "subscriber never called"
        
    def test_call_depth(self):
         
        def recurrent_method(_):
            self.event_broker.publish(DummyEvent())
            
        self.event_broker.register(DummyEvent, recurrent_method)
        
        with pytest.raises(RuntimeError, match=f"Event recursive depth of {self.event_broker.call_depth_max} has been exceeded."):
            self.event_broker.publish(DummyEvent())

    def test_call_depth_does_not_leak_when_subscriber_raises(self):
        """A subscriber exception must not permanently inflate call_depth (it's a Singleton with no reset)."""

        def raises(_):
            raise ValueError("boom")

        self.event_broker.register(DummyEvent, raises)

        with pytest.raises(ValueError, match="boom"):
            self.event_broker.publish(DummyEvent())

        assert self.event_broker.call_depth == 0

        # A subsequent, unrelated publish must still work rather than immediately
        # hitting the recursion-depth ceiling from the earlier leaked call_depth.
        method_called = False

        def dummy_method(_):
            nonlocal method_called
            method_called = True

        self.event_broker.register(SmartyEvent, dummy_method)
        self.event_broker.publish(SmartyEvent())
        assert method_called

    def test_multiple_subscribers(self):
        e = DummyEvent()
        
        method1_called = False
        def dummy_method_1(_):
            nonlocal method1_called
            method1_called = True
            
        method2_called = False
        def dummy_method_2(_):
            nonlocal method2_called
            method2_called = True
            
        method3_called = False
        def smarty_method(event :SmartyEvent):
            nonlocal method3_called
            method3_called = True
            assert type(event) is SmartyEvent, "wrong event type forwarded to subscriber"
            
        self.event_broker.register(DummyEvent, dummy_method_1)
        self.event_broker.register(DummyEvent, dummy_method_2)
        self.event_broker.register(SmartyEvent, smarty_method)
        
        self.event_broker.publish(e)
        
        assert not method3_called, "Incorrectly called different event"
        
        assert method1_called, "Failed to call subscriber 1"
        assert method2_called, "Failed to call subscriber 2"
        
        self.event_broker.publish(SmartyEvent())
        assert method3_called, "Failed to call subscriber 3"