"""Decision policy injection; Gaussian control always goes through Session tools."""
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from typing import Protocol
from .models import finite, integer

SYSTEM_PROMPT = '''You are a computational chemistry research agent using Gaussian.
Choose exactly ONE next tool action based on the observation, current geometry,
bond changes, previous failures, job diagnostics and remaining budgets. Your
choices need not follow a fixed pipeline. Preserve stable mapped atom IDs.
Return JSON: {"tool": "...", "params": {...}, "reason": "chemical and numerical justification"}.
Use only the tools listed in the observation. Never return shell commands, code,
raw Gaussian routes or coordinates masquerading as validated structures.
Geometry tools: set_distance(atoms=[anchor,target], value=Angstrom, moving=[IDs]);
set_angle(atoms=[a,b,c], value=degrees, moving=[IDs including c]);
set_dihedral(atoms=[a,b,c,d], value=degrees, moving=[IDs including d]);
translate_fragment(moving=[IDs], vector=[dx,dy,dz]);
rotate_fragment(moving=[IDs], origin=[x,y,z], axis=[x,y,z], degrees=number).
Moving fragments must exclude reference anchor atoms; edits are rigid.
run_gaussian kind is opt, frozen_opt, scan, ts, qst2, qst3, freq,
irc_forward, irc_reverse, endpoint_forward or endpoint_reverse.
For frozen_opt/scan constraints use {"type":"B|A|D", "atoms":[IDs],
"operation":"F|S", "steps":integer, "step_size":number}; steps/step_size only for S.
Frequency/IRC need the TS output structure_id and parent_job; endpoint optimization
needs the corresponding IRC output structure_id and parent_job. All child jobs
must use the same method, basis, solvent and integration grid as their parent.
Use use_structure(structure_id) to branch or undo_geometry_edit() to revert.
NoEigenTest is not available. Use request_human_review(question) for ambiguity.
A normal termination or one imaginary frequency NEVER proves the target TS.
Only verify_ts(ts_job_id optional) may return VALIDATED_TS. It requires
unconstrained TS convergence, significant single negative curvature, target
mode analysis, two complete IRC paths and optimized matching mapped endpoints.
Synthetic/offline evidence cannot be validated. Do not claim success yourself.
Stop by pause() or request_human_review(), preserving uncertainty and rationale.'''


class DecisionPolicy(Protocol):
    def decide(self, observation):
        """Return one action dict or None when no further decision is available."""


class ReplayPolicy:
    """Explicit human action plans for reproducibility; not an autonomous LLM."""
    def __init__(self, actions):
        if not isinstance(actions, list):
            raise ValueError('replay plan must be a JSON list')
        self.actions = iter(actions)

    def decide(self, observation):
        return next(self.actions, None)


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError('LLM endpoint redirects are disabled')


class HTTPPolicy:
    """User-configured chat-completions endpoint. API keys stay in environment."""
    def __init__(self, endpoint, model, api_key_env='TS_LLM_API_KEY', timeout_seconds=45., max_tokens=1500):
        url = urllib.parse.urlsplit(endpoint)
        localhost = url.hostname in ('127.0.0.1', 'localhost', '::1')
        if not url.hostname or (url.scheme != 'https' and not (url.scheme == 'http' and localhost)) or url.username or url.password or url.query or url.fragment:
            raise ValueError('LLM endpoint must use HTTPS (HTTP only for loopback), without embedded credentials or query')
        if not isinstance(model, str) or not model.strip():
            raise ValueError('LLM model is required')
        self.timeout = finite(timeout_seconds, 'LLM timeout')
        if not 0 < self.timeout <= 60:
            raise ValueError('LLM timeout must be in (0,60] seconds')
        integer(max_tokens, 'LLM max_tokens', 1)
        if api_key_env is not None and (not isinstance(api_key_env, str) or not api_key_env):
            raise ValueError('api_key_env must be an environment variable name or null')
        key = os.environ.get(api_key_env) if api_key_env else None
        if api_key_env and not key:
            raise ValueError('missing API key environment variable: %s' % api_key_env)
        self.endpoint, self.model, self.key, self.max_tokens = endpoint, model, key, max_tokens

    def decide(self, observation):
        payload = {'model': self.model, 'messages': [{'role': 'system', 'content': SYSTEM_PROMPT},
                                                   {'role': 'user', 'content': json.dumps(observation, ensure_ascii=False, allow_nan=False)}],
                   'response_format': {'type': 'json_object'}, 'max_tokens': self.max_tokens}
        headers = {'Content-Type': 'application/json'}
        if self.key:
            headers['Authorization'] = 'Bearer ' + self.key
        req = urllib.request.Request(self.endpoint, data=json.dumps(payload).encode(), headers=headers, method='POST')
        try:
            with urllib.request.build_opener(NoRedirect()).open(req, timeout=self.timeout) as response:
                raw = response.read(1048577)
                if len(raw) > 1048576:
                    raise ValueError('LLM response exceeds 1 MiB')
        except urllib.error.HTTPError as exc:
            raise ValueError('LLM HTTP error %d; response body omitted' % exc.code) from None
        except (urllib.error.URLError, TimeoutError, OSError):
            raise ValueError('LLM endpoint unavailable or timed out') from None
        try:
            data = json.loads(raw)
            content = data['choices'][0]['message']['content']
            action = json.loads(content)
            if not isinstance(action, dict):
                raise ValueError('action must be object')
            return action
        except (ValueError, TypeError, KeyError, IndexError):
            raise ValueError('LLM response must contain a single JSON action in choices[0].message.content') from None


class Agent:
    def __init__(self, session, policy, on_step=None):
        self.session, self.policy, self.on_step = session, policy, on_step

    def run(self, max_steps=20):
        integer(max_steps, 'agent max_steps', 1)
        results, stop = [], 'step_limit'
        for _ in range(max_steps):
            observation = self.session.observe()
            if observation['paused']:
                stop = 'paused'
                break
            if observation['usage']['actions'] >= observation['limits']['max_actions']:
                stop = 'action_budget'
                break
            try:
                action = self.policy.decide(observation)
            except (ValueError, OSError) as exc:
                self.session.execute({'tool': 'request_human_review', 'params': {'question': str(exc)},
                                      'reason': 'decision provider failed; human intervention required', 'actor': 'agent'})
                stop = 'policy_error'
                break
            if action is None:
                stop = 'policy_exhausted'
                break
            if isinstance(action, dict):
                action = dict(action, actor='agent')  # A model cannot impersonate a human.
            result = self.session.execute(action)
            results.append(result)
            if self.on_step:
                self.on_step(action, result)
            if result['ok'] and isinstance(result.get('result'), dict):
                if result['result'].get('label') == 'VALIDATED_TS':
                    stop = 'validated_ts'
                    break
                if self.session.state['paused']:
                    stop = 'human_review' if action.get('tool') == 'request_human_review' else 'paused'
                    break
        report = self.session.evaluate()
        return {'stop_reason': stop, 'steps': len(results), 'actions': results, 'report': report}
