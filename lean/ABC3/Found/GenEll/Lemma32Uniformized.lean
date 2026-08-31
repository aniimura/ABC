/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma32StableLine

/-!
# [GenEll] Lemma 3.2, (i) —— **同変な一意化を入力として**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★★逸脱が 1 つ消えた

`Lemma32StableLine.lean` の `lemma_3_2_i` は `.src` にこう書いていた:

> `Lemma 3.2, (i)（Kˣ/q^ℤ の側——**`G_K` 同変な一意化は含まない**）`

★★**その一意化は 2026-08-27 に取れた**（`Found/GaloisRep/` の 8 段、
`tatePhi_map`・`tatePhiAddEquivAll`）。
★本ファイルはそれを**入力として受ける形**に書き直し、
`M_l(E)` の側から `Kˣ/q^ℤ` の側へ渡す橋を架ける。

## ★★★★★★橋の中身

`M_l(E)` の `G_K`-安定な直線 `⟨P⟩` から `lemma_3_2_i` の仮説を作る:

1. `c ≝ Φ⁻¹(P)` は `Lˣ/Q^ℤ` の `l`-捩れ。
2. `c` を `x ∈ Lˣ` へ持ち上げると `x^l = Q^k`（`exists_lift_pow`）。
3. 安定性 `σP = c′·P` に `Φ⁻¹` を当てると `mk(σx) = mk(x^{c′})`、
   したがって **`σx = x^{c′}·Q^n`** ——これが `lemma_3_2_i` の `hstab` である。

★★機構は `Φ` が**同変かつ単射**であることだけである。

## ★残っている段（明示）

★`Lemma 3.2` 全体には `(ii)`（`E/μ_l` が母数 `q^l` の Tate 曲線）も要る。
★★`Lemma32QuotMu.lean` が**一意化の側で**それを持っているが、
`E/μ_l` を**スキームとして**作る段（有限平坦群スキームによる商）は
mathlib にも ABC3 にも無い。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep

/-! ## ★捩れ元の持ち上げ -/

/-- ★★**`l`-捩れの類は `x^l = Q^k` の形に持ち上がる**。

★`QuotientGroup.mk` は全射なので代表 `x` が取れ、
`mk (x^l) = 1` から `x^l ∈ zpowers Q` が出る。 -/
theorem exists_lift_pow {L : Type} [Field L] (Q : Lˣ) (l : ℕ)
    (c : Lˣ ⧸ Subgroup.zpowers Q) (hc : c ^ l = 1) :
    ∃ (x : Lˣ) (k : ℤ), (QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers Q) = c ∧ x ^ l = Q ^ k := by
  obtain ⟨x, rfl⟩ := QuotientGroup.mk_surjective c
  have hmk : (QuotientGroup.mk (x ^ l) : Lˣ ⧸ Subgroup.zpowers Q)
      = (QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers Q) ^ l := rfl
  have hmem : x ^ l ∈ Subgroup.zpowers Q := by
    refine (QuotientGroup.eq_one_iff (x ^ l)).1 ?_
    rw [hmk]
    exact hc
  obtain ⟨k, hk⟩ := hmem
  exact ⟨x, k, rfl, hk.symm⟩

/-! ## ★★★★★★★★★★★同変な一意化からの `Lemma 3.2, (i)` -/

/-- ★★★★★★★★★★★**[GenEll] Lemma 3.2, (i) —— 同変な一意化を入力として**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

`M_l(E)` の `G_K`-安定な直線 `⟨P⟩` が `𝔽_l(1)` でない（`¬ l ∣ k`）なら
**`l ∣ v_K(q_E)`** である。

★★入力は「**同変な一意化** `Φ : Lˣ/Q^ℤ ≃+ E`」だけである
——それは `Found/GaloisRep/` の 8 段（2026-08-27）で構成した:
`tatePhiAddEquivAll`（全単射）＋ `tatePhi_map`（同変性）。

★★★これで `Lemma32StableLine.lean` の逸脱
「`G_K` 同変な一意化は含まない」が**消える**。

★★★★機構は `Φ` が同変かつ単射であることだけである。 -/
theorem lemma_3_2_i_of_uniformization {K L : Type} [Field K] [Field L] [Algebra K L] [IsGalois K L]
    {l : ℕ} (hl : l.Prime) (q : Kˣ) (v : Kˣ →* Multiplicative ℤ)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0)
    {E : Type} [AddCommGroup E]
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q)) ≃+ E)
    (act : (L ≃ₐ[K] L) → E →+ E)
    (hequiv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ),
      act σ (Φ (Additive.ofMul (QuotientGroup.mk u)))
        = Φ (Additive.ofMul (QuotientGroup.mk (Units.map (σ : L →* L) u))))
    (P : E) (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ, act σ P = c • P)
    (x : Lˣ) (k : ℤ)
    (hxP : Φ (Additive.ofMul (QuotientGroup.mk x)) = P)
    (hxk : x ^ l = (Units.map (algebraMap K L : K →* L) q) ^ k)
    (hk : ¬ ((l : ℤ) ∣ k)) :
    (l : ℤ) ∣ vAdd v q := by
  refine lemma_3_2_i hl q v hqinf x k hxk hk (fun σ => ?_)
  obtain ⟨c, hc⟩ := hstab σ
  have h1 : Φ (Additive.ofMul (QuotientGroup.mk (Units.map (σ : L →* L) x)))
      = Φ (Additive.ofMul
        ((QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q))
          ^ c)) := by
    rw [← hequiv σ x, hxP, hc, ← hxP, ← map_zsmul]
    rfl
  have h3 : (QuotientGroup.mk (Units.map (σ : L →* L) x)
        : Lˣ ⧸ Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q))
      = QuotientGroup.mk (x ^ c) := Additive.ofMul.injective (Φ.injective h1)
  have h4 : (Units.map (σ : L →* L) x) * (x ^ c)⁻¹
      ∈ Subgroup.zpowers (Units.map (algebraMap K L : K →* L) q) := by
    rw [← QuotientGroup.eq_one_iff, QuotientGroup.mk_mul, h3, QuotientGroup.mk_inv,
      mul_inv_cancel]
  obtain ⟨n, hn⟩ := h4
  refine ⟨c, n, ?_⟩
  have hn' : (Units.map (algebraMap K L : K →* L) q) ^ n
      = (Units.map (σ : L →* L) x) * (x ^ c)⁻¹ := hn
  rw [hn', ← mul_assoc, mul_comm (x ^ c) _, mul_assoc, mul_inv_cancel, mul_one]

/-! ### ★出典の紐付け(`.src`)

★★★**逸脱「`G_K` 同変な一意化は含まない」が消えた。**
★ただし `Lemma 3.2` 全体には `(ii)` も要り、
そこには `E/μ_l` を**スキームとして**作る段が残っている。 -/

def exists_lift_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(l-捩れの類を x^l = Q^k の形へ持ち上げる)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2_i_of_uniformization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(同変な一意化を入力として。(ii) は含まない)",
    sectionId := "genell-lemma-3-2" }

def lemma_3_2_i_of_uniformization.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_2_i(Kˣ/q^ℤ の側の二者択一)"
      (.inProject "ABC3" "ABC3.Found.GenEll.lemma_3_2_i") 15,
    .citation "[ABC3]" "tatePhi_map(Φ の G_K-同変性——2026-08-27 の 8 段)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_map") 15,
    .citation "[ABC3]" "tatePhiAddEquivAll(Φ が全単射であること)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhiAddEquivAll") 15,
    .otherPaper "[FC]" "Degenerations of Abelian Varieties, Chapter III, Corollary 7.3" 15,
    .implicitStep
      ("★★★★★**逸脱「G_K 同変な一意化は含まない」が消えた**" ++
       "——mathlib-gap の galois-equivariant-tate-uniformization は" ++
       "2026-08-27 に 10 段中 9 段まで閉じ、本ファイルが橋を架けた") 15,
    .implicitStep
      ("★★残る段: Lemma 3.2 全体には (ii)(E/μ_l が母数 q^l の Tate 曲線)も要る。" ++
       "★Lemma32QuotMu.lean が**一意化の側で**それを持っているが、" ++
       "E/μ_l を**スキームとして**作る段(有限平坦群スキームによる商)は" ++
       "mathlib にも ABC3 にも無い") 15 ]

end ABC3.Found.GenEll
