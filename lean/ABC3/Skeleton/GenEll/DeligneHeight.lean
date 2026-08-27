/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.FiniteFromNorthcott

/-!
# [DelSB616] Exemple 1.4 —— **`[Silv2]` の代わりになる原典**（`Skeleton`）

原典: P. Deligne, *Preuve des conjectures de Tate et de Shafarevitch (d'après G. Faltings)*,
Séminaire Bourbaki 616 (1983/84)、物理 p.5。

原文 (DelSB616 p.5):
> Exemple 1.4 Le cas k = Q , g = 1 .- Soient E une courbe elliptique sur Q et w

## ★★★★★★★★なぜこれを載せるのか —— 2026-08-27 の訂正

`[GenEll] Proposition 3.4` の不等式 `deg ≾ ht_∞ ≾ 12(1+ε)·ht_Falt ≾ (1+ε)·ht_∞` は
原文が `[Silv2]`（Silverman, *Advanced Topics*）Prop 2.1 に帰している。
★**`[Silv2]` は `0_Source` に無い**ので「原典が無い」と台帳に書いていた。

★★★**それは誤りだった。** `[DelSB616]` §1 が同じ内容を扱っており、
**`0_Source` にあり `papers.json` に登記済み**である。
★Deligne 本人が「C'est la partie la plus technique, et la plus longue, de la preuve」
と書いている段が、まさにこれである。

★★★★したがって CLAUDE.md の進め方（スケルトンで依存グラフを作成後、葉から形式化）が
**そのまま適用できる**。「壁」ではなく「まだグラフに載せていない道」である。

## ★★★★★★目視で確認した 6 段（260 dpi、2026-08-27）

`H(E) ≝ exp h(E)`（Faltings 高さ）、`H*(E) ≝ sup(|c₄|^{1/4}, |c₆|^{1/6})`（素朴な高さ）として:

| 段 | 内容 | 見込み |
|---|---|---|
| 1 | `(E,ω) ≅ (ℂ/Λ, dz)`、`Λ = ℤω₁ ⊕ ℤω₂`、`τ = ω₂/ω₁`、`q = exp(2πiτ)` | 複素解析（並行セッションが `LatticeCurve.lean` 等で構築中） |
| 2 | `H(E) = (1/|ω₁|)·|Im(τ)/π|^{-1/2}`（`|ω₁|` 最小の基底で） | 同上 |
| 3 | `H*(E) ≤ C₀·H(E)·Im(τ)^{1/2}` | 同上 |
| 4 | `j = c₄³/Δ = q⁻¹ + 744 + …`、`Im(τ) ∼ sup(1, log|j|)/2π`（`Im(τ) > a`） | ★`JFunction.lean` に `tendsto_norm_jFun_atImInfty` あり |
| 5 | `Δ` は整数なので `|j| ≤ H*(E)^{12}`、ゆえに `H*(E) ≤ C₁·H(E)·sup(1, log H*(E))^{1/2}` | 初等 |
| 6 | ★**`H*(E) ≤ C₂·H(E)·sup(1, log H(E))^{1/2}`** | ★**初等**（`Lemma 4.2` と同型の評価） |

★★段 5 → 6 は原文 `[GenEll] Lemma 4.2`（An Elementary Estimate）と**同じ形**であり、
`Found/GenEll/Elementary.lean` に既にある。

## ★★★★★★★★これが依存グラフのどこに入るか

    [DelSB616] Exemple 1.4  ⟶  [GenEll] Proposition 3.4 の不等式
                            ⟶  Lemma 3.5 / Lemma 3.7 / Theorem 3.8 / Cor 4.3 / Cor 4.4

★`Proposition 3.4` の**有限性の主張**のほうは既に取れている
（`Found/GenEll/FiniteFromNorthcott.lean` の `finite_of_le_of_northcott`）。
★★**残るのは不等式の側だけ**であり、その原典が本ファイルの節点である。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Found.GenEll

/-- ★★★★★★★★**[DelSB616] Exemple 1.4 の到達点** ——
素朴な高さは Faltings 高さで（対数の補正つきで）抑えられる。

原文 (DelSB616 p.5):
> Exemple 1.4 Le cas k = Q , g = 1 .- Soient E une courbe elliptique sur Q et w

★対数を取った形で述べる: `log H* ≤ log H + (1/2)·log(sup(1, log H)) + C`。

★★★**これが `[GenEll] Proposition 3.4` の `ht_∞ ≾ a·ht_Falt` に当たる**
——`Found/GenEll/FiniteFromNorthcott.lean` の `finite_of_le_of_northcott` の
第 1 の入力である。

★★段 6 の形は `[GenEll] Lemma 4.2`（An Elementary Estimate）と同型で、
`Found/GenEll/Elementary.lean` に既にある。 -/
theorem deligne_height_comparison
    {P : Type*} (logHstar logH : P → ℝ)
    (hH : ∀ p, 0 ≤ logH p)
    -- ★段 3（複素解析）: `H* ≤ C₀·H·Im(τ)^{1/2}`
    (imTau : P → ℝ) (hIm : ∀ p, 0 < imTau p) (C₀ : ℝ)
    (h3 : ∀ p, logHstar p ≤ logH p + (1/2) * Real.log (imTau p) + C₀)
    -- ★段 4（`j` の漸近）＋段 5（`Δ` が整数）: `Im(τ) ≤ C₁·sup(1, log H*)`
    (C₁ : ℝ) (hC₁ : 1 ≤ C₁)
    (h45 : ∀ p, imTau p ≤ C₁ * max 1 (logHstar p)) :
    ∃ C₂ : ℝ, ∀ p, logHstar p ≤ logH p + (1/2) * Real.log (max 1 (logHstar p)) + C₂ := by
  refine ⟨C₀ + (1/2) * Real.log C₁, fun p => ?_⟩
  have hpos : (0:ℝ) < max 1 (logHstar p) := lt_of_lt_of_le zero_lt_one (le_max_left _ _)
  have hle : Real.log (imTau p) ≤ Real.log C₁ + Real.log (max 1 (logHstar p)) := by
    calc Real.log (imTau p) ≤ Real.log (C₁ * max 1 (logHstar p)) :=
          Real.log_le_log (hIm p) (h45 p)
      _ = Real.log C₁ + Real.log (max 1 (logHstar p)) :=
          Real.log_mul (by linarith) (ne_of_gt hpos)
  have := h3 p
  linarith

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def deligne_height_comparison.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 5,
    item := "Exemple 1.4(素朴な高さは Faltings 高さで抑えられる——段 3・段 4-5 は仮定に置く)",
    sectionId := "delsb616-ex-1-4" }

/-- ★原文 p.5 を 260 dpi で目視して数えた（2026-08-27）。 -/
def deligne_height_comparison.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★段 1-2(複素解析): (E,ω) ≅ (ℂ/Λ, dz)、|ω₁| 最小の基底で " ++
       "H(E) = (1/|ω₁|)·|Im(τ)/π|^{-1/2}。並行セッションが LatticeCurve.lean / " ++
       "Covolume.lean 等で構築中") 5,
    .implicitStep
      ("★段 3(複素解析): H*(E) ≤ C₀·H(E)·Im(τ)^{1/2}。本 statement は仮定 h3 に置いた") 5,
    .implicitStep
      ("★段 4: j = c₄³/Δ = q⁻¹ + 744 + … と Im(τ) ∼ sup(1,log|j|)/2π(Im(τ) > a)。" ++
       "★Found/GenEll/JFunction.lean の tendsto_norm_jFun_atImInfty が入口") 5,
    .implicitStep
      ("★段 5: Δ は整数なので |j| ≤ H*(E)^{12}。初等。段 4 と合わせて h45 に置いた") 5,
    .implicitStep
      ("★★段 6: H*(E) ≤ C₂·H(E)·sup(1, log H(E))^{1/2}。" ++
       "★[GenEll] Lemma 4.2(An Elementary Estimate)と同型の評価で、" ++
       "Found/GenEll/Elementary.lean に既にある。本 statement は log H* の側で止めており、" ++
       "log H の側へ移す段が残っている") 5,
    .otherPaper "[GenEll]"
      ("Proposition 3.4 の不等式 deg ≾ ht_∞ ≾ 12(1+ε)·ht_Falt ≾ (1+ε)·ht_∞ —— " ++
       "★本節点がその原典である。原文は [Silv2] Prop 2.1 に帰すが [Silv2] は 0_Source に" ++
       "無く、[DelSB616] が同じ内容を扱う") 17,
    .implicitStep
      ("★★★★★2026-08-27 の訂正: 「mathlib に無い」と「原典が無い」を混同していた。" ++
       "[Silv2] は手元に無いが [DelSB616] / [Elkik] / [RayHt] が手元にあり登記済みである。" ++
       "CLAUDE.md の進め方(スケルトンで依存グラフを作成後、葉から形式化)がそのまま適用できる") 5 ]

end ABC3.Skeleton.GenEll
