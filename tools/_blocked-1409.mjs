import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const key = '[GenEll] Lemma 3.5(Vélu の商の半安定性——判別式の恒等式と p ∣ l の形式群)';
const e = j.blocked.find(b => b.key === key);
if (!e) { console.log('not found'); process.exit(1); }
e.nonMuKernelRoute2026_09_02 =
  '★★★★★★★★**`l ∣ v_P(q)` の場合の道が見えた（2026-09-02、第 1409 の測定）**'
  + '——`hcop` は**仮定から外せる**見込みである。'
  + '☆Tate 母数 `Q` に対し、位数 `l` の点の代表 `u ∈ Kˣ` は `u^l = Q^m` を満たす。'
  + '★**二者択一**:（i）`l ∣ m` なら `u·Q^{−m/l}` に取り替えて `u^l = 1`'
  + '（＝`μ_l` 型、第 1404-1406 の道。`hcop` はこれを強制するために使っていた）。'
  + '（ii）`l ∤ m` なら `vAdd(u) = m·vAdd(Q)/l` は `vAdd(Q)ℤ` に入らないので、'
  + '代表を `0 < vAdd(u) < vAdd(Q)` に取れば **`u ∈ 𝔪`**。'
  + '☆(ii) では核の非零点すべてが `u^k`（`1 ≤ k ≤ l−1`、代表の付値は `(bk mod l)·vAdd(Q)/l > 0`）'
  + 'で `𝔪` に入るので、**`x` 座標がすべて `𝔪` に入る**'
  + '（`tateXpair a b` は `a, b ∈ 𝔪` なら `a/(1−a)² + O(q) ≡ 0`）。'
  + '★したがって `veluV2 = 6x² + b₂x + b₄ ≡ b₄ = −10s₃(q) ≡ 0 (mod 𝔪)`、'
  + 'すなわち `v = Σ veluV2 ∈ 𝔪` で '
  + '**`c₄(veluCurve) = c₄(E_q) + 240v ≡ c₄(E_q)` は単元**——半安定性が出る。'
  + '☆(i) の `μ_l` 型では逆に `c₄′ ≡ l⁴` が単元（第 1388 の `h4`）。'
  + '★★★どちらの場合も `c₄′` が単元なので、二者択一を作れば'
  + '**悪い素点は `hcop` なしで完全に閉じる**。'
  + '☆見積もり: 二者択一（1-2）＋ `tateXpair ∈ 𝔪` の評価（1-2）＋ Vélu 和と `c₄`（1-2）'
  + '＋ 組み立て（1）＝ **4-7 ブロック**。'
  + '★これは `hcopDoesNotDescend2026_09_02` の問題（底変換で `l ∤ jExp` が壊れる）を'
  + '**根本から回避する**——そもそも仮定しなくてよくなる。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded blocked-leaves');
