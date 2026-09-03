/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TatePhi
import ABC3.Found.GaloisRep.TateClassPt

/-!
# Galois (G6) 第 884 ブロック —— **★★★★★★★★★★`μ_l` の点は `tatePtPair ζⁱ (qζ⁻ⁱ)`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★これは何か

`Lemma 3.5` の葉 1 に残っているのは、「`H` が各悪い素点で `μ_l` に
対応する」という仮説から `hquot`（Vélu の商の `j`）を導く**配管**である。

★その第一歩が本ファイルである:

    一意化 `Φ : Kˣ/q^ℤ ≃ E_q(K)` による `ζⁱ` の行き先は
    **`tatePtPair ζⁱ (q·ζ⁻ⁱ)`**、すなわち `tateXpair (ζⁱ) (q(ζⁱ)^{l-1})` である。

☆これは `c4_velu_tate`・`c6_velu_tate`（第 853・867）が使っている座標そのものであり、
両者を繋ぐのが本ブロックの役割である。

★機構は一本道である。`normRep`（基本領域の代表元）は
`0 ≤ v(u) < v(q)` で決まるが、1 の冪根は `v = 0` なので**自分自身が代表元**である。
したがって `tateAOf S [ζⁱ] = ζⁱ` であり、`a·w = q` から `tateWOf S [ζⁱ] = q·(ζⁱ)^{l-1}` である。

| 定理 | 内容 |
|---|---|
| `vAdd_eq_zero_of_pow_eq_one` | ★1 の冪根の付値は `0` |
| `normRep_of_pow_eq_one` | ★★1 の冪根の類の代表元はそれ自身 |
| `tateAOf_of_pow_eq_one` | ★★★★`tateAOf S [ζⁱ] = ζⁱ` |
| `tateWOf_of_pow_eq_one` | ★★★★`tateWOf S [ζⁱ] = q·(ζⁱ)^{l-1}` |
| `tatePhi_of_pow_eq_one` | ★★★★★★★★★★**`Φ(ζⁱ) = tatePtPair ζⁱ (q(ζⁱ)^{l-1})`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★1 の冪根の付値 -/

variable {K : Type} [Field K]

/-- ★**1 の冪根の付値は `0`**——`l·v(u) = v(u^l) = v(1) = 0` だから。 -/
theorem vAdd_eq_zero_of_pow_eq_one (v : Kˣ →* Multiplicative ℤ) {l : ℕ} (hl : 0 < l)
    (u : Kˣ) (h : u ^ l = 1) : vAdd v u = 0 := by
  have h1 : (l : ℤ) * vAdd v u = 0 := by
    rw [← vAdd_zpow v u (l : ℤ), zpow_natCast, h, vAdd_one]
  have hl0 : (l : ℤ) ≠ 0 := by exact_mod_cast hl.ne'
  exact (mul_eq_zero.1 h1).resolve_left hl0

/-- ★★**1 の冪根の類の正規化代表元はそれ自身**。 -/
theorem normRep_of_pow_eq_one (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {l : ℕ} (hl : 0 < l) (u : Kˣ) (h : u ^ l = 1) :
    normRep v Q hQ (QuotientGroup.mk u) = u := by
  refine normRep_eq_self v Q hQ u ?_ ?_
  · rw [vAdd_eq_zero_of_pow_eq_one v hl u h]
  · rw [vAdd_eq_zero_of_pow_eq_one v hl u h]
    exact hQ

/-! ## ★★★★正規化した対 -/

section Pair

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [Algebra R K]

/-- ★★★★**`tateAOf S [ζⁱ] = ζⁱ`**——付値が `0` なので代表元はそれ自身。 -/
theorem tateAOf_of_pow_eq_one (S : TateSetup R I K) {l : ℕ} (hl : 0 < l)
    (a : R) (u : Kˣ) (hau : algebraMap R K a = (u : K)) (hul : u ^ l = 1) :
    tateAOf S (QuotientGroup.mk u) = a := by
  refine S.hinj ?_
  rw [(tateAOf_spec S (QuotientGroup.mk u)).1,
    normRep_of_pow_eq_one S.v S.Q S.hQ hl u hul, hau]

/-- ★`a` の側も `l` 乗すると `1` になる。 -/
theorem pow_eq_one_of_map (S : TateSetup R I K) {l : ℕ}
    (a : R) (u : Kˣ) (hau : algebraMap R K a = (u : K)) (hul : u ^ l = 1) :
    a ^ l = 1 := by
  refine S.hinj ?_
  rw [map_pow, hau, map_one]
  exact congrArg (fun x : Kˣ => (x : K)) hul

/-- ★★★★**`tateWOf S [ζⁱ] = q·(ζⁱ)^{l-1}`**——`a·w = q` と `a^l = 1` から。 -/
theorem tateWOf_of_pow_eq_one (S : TateSetup R I K) {l : ℕ} (hl : 0 < l)
    (a : R) (u : Kˣ) (hau : algebraMap R K a = (u : K)) (hul : u ^ l = 1) :
    tateWOf S (QuotientGroup.mk u) = S.q * a ^ (l - 1) := by
  have hal : a ^ l = 1 := pow_eq_one_of_map S a u hau hul
  have hsplit : a * a ^ (l - 1) = 1 := by
    have hl' : l - 1 + 1 = l := Nat.succ_pred_eq_of_pos hl
    calc a * a ^ (l - 1) = a ^ (l - 1 + 1) := by ring
      _ = a ^ l := by rw [hl']
      _ = 1 := hal
  have hunit : IsUnit a :=
    ⟨⟨a, a ^ (l - 1), hsplit, by rw [mul_comm]; exact hsplit⟩, rfl⟩
  have hmul := tateAOf_mul_tateWOf S (QuotientGroup.mk u)
  rw [tateAOf_of_pow_eq_one S hl a u hau hul] at hmul
  refine hunit.mul_left_cancel ?_
  rw [hmul]
  linear_combination -S.q * hsplit

end Pair

/-! ## ★★★★★★★★★★`Φ` の行き先 -/

section Phi

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [DecidableEq K] [Algebra R K]

/-- ★`tatePtPair` は値だけで決まる（証明の引数は無関係）。

☆`rw` だと motive が壊れるので、**`subst` で潰せる形**にしておく。 -/
theorem tatePtPair_congr {a w a' w' q : R} (hq : q ∈ I)
    {haw : a * w = q} {hw : IsUnit (1 - w)} {hne : algebraMap R K (1 - a) ≠ 0}
    {haw' : a' * w' = q} {hw' : IsUnit (1 - w')} {hne' : algebraMap R K (1 - a') ≠ 0}
    {hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0}
    (ha : a = a') (hww : w = w') :
    tatePtPair a w q hq haw hw hne hΔ = tatePtPair a' w' q hq haw' hw' hne' hΔ := by
  subst ha
  subst hww
  rfl

/-- ★★★★★★★★★★**`Φ(ζⁱ) = tatePtPair ζⁱ (q(ζⁱ)^{l-1})`**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★これが `μ_l ⊂ Kˣ/q^ℤ` と `c4_velu_tate`・`c6_velu_tate` の座標を繋ぐ一本である。 -/
theorem tatePhi_of_pow_eq_one (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : 0 < l) (a : R) (u : Kˣ)
    (hau : algebraMap R K a = (u : K)) (hul : u ^ l = 1)
    (hc : (QuotientGroup.mk u : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1)
    (haw : a * (S.q * a ^ (l - 1)) = S.q)
    (hw : IsUnit (1 - S.q * a ^ (l - 1)))
    (hne : algebraMap R K (1 - a) ≠ 0) :
    tatePhi S hΔ (QuotientGroup.mk u)
      = tatePtPair a (S.q * a ^ (l - 1)) S.q S.hq haw hw hne hΔ := by
  rw [tatePhi_eq S hΔ hc]
  exact tatePtPair_congr S.hq
    (tateAOf_of_pow_eq_one S hl a u hau hul)
    (tateWOf_of_pow_eq_one S hl a u hau hul)

end Phi

/-! ## ★★★★★★類の上でも `μ_l` は潰れない -/

/-- ★★★**1 の冪根の類が単位ならそれ自身が `1`**。

☆`u = Q^n` と書けるなら `0 = v(u) = n·v(Q)` で、`v(Q) > 0` なので `n = 0`。 -/
theorem eq_one_of_mk_eq_one_of_pow_eq_one (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ)
    (hQ : 0 < vAdd v Q) {l : ℕ} (hl : 0 < l) (u : Kˣ) (hul : u ^ l = 1)
    (h : (QuotientGroup.mk u : Kˣ ⧸ Subgroup.zpowers Q) = 1) : u = 1 := by
  obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff u).1 h
  have hn' : Q ^ n = u := hn
  have hv : vAdd v u = 0 := vAdd_eq_zero_of_pow_eq_one v hl u hul
  rw [← hn', vAdd_zpow] at hv
  have hn0 : n = 0 := by
    rcases mul_eq_zero.1 hv with h1 | h1
    · exact h1
    · omega
  rw [← hn', hn0, zpow_zero]

/-- ★★★★★★**`μ_l` の元は `Kˣ/Q^ℤ` の中でも区別される**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆`v(q) > 0` なので `q^ℤ` は 1 の冪根を 1 しか含まない。 -/
theorem mk_pow_injOn (v : Kˣ →* Multiplicative ℤ) (Q : Kˣ) (hQ : 0 < vAdd v Q)
    {l : ℕ} (hl : 0 < l) (ζ : Kˣ) (hζl : ζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → ζ ^ n ≠ 1)
    {i j : ℕ} (hi : i < l) (hj : j < l)
    (h : (QuotientGroup.mk (ζ ^ i) : Kˣ ⧸ Subgroup.zpowers Q)
       = QuotientGroup.mk (ζ ^ j)) : i = j := by
  -- ★差の側を見る。`j ≤ i` の場合を補顓にして 2 回使う
  have key : ∀ {i j : ℕ}, i < l → j ≤ i →
      (QuotientGroup.mk (ζ ^ i) : Kˣ ⧸ Subgroup.zpowers Q)
        = QuotientGroup.mk (ζ ^ j) → i = j := by
    intro i j hi hji h
    have hsplit : ζ ^ i = ζ ^ (i - j) * ζ ^ j := by
      rw [← pow_add]
      congr 1
      omega
    have hone : (QuotientGroup.mk (ζ ^ (i - j)) : Kˣ ⧸ Subgroup.zpowers Q) = 1 := by
      have h3 : (QuotientGroup.mk (ζ ^ (i - j)) : Kˣ ⧸ Subgroup.zpowers Q)
          = QuotientGroup.mk (ζ ^ i) * (QuotientGroup.mk (ζ ^ j))⁻¹ := by
        rw [← QuotientGroup.mk_inv, ← QuotientGroup.mk_mul]
        congr 1
        rw [hsplit]
        group
      rw [h3, h, mul_inv_cancel]
    have hpow : (ζ ^ (i - j)) ^ l = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, hζl, one_pow]
    have hz := eq_one_of_mk_eq_one_of_pow_eq_one v Q hQ hl (ζ ^ (i - j)) hpow hone
    by_contra hne
    exact hord (i - j) (by omega) (by omega) hz
  rcases le_total j i with hji | hij
  · exact key hi hji h
  · exact (key hj hij h.symm).symm

/-! ## ★★★★★★★★`l ∣ k` なら類は 1 の冪根の類である -/

/-- ★★★★★★★★**`x^l = Q^k` で `l ∣ k` なら、`[x]` は `μ_l` の元の類である**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★これが `Lemma 3.2, (i)` の**結論側**である——
`lemma_3_2_i_tate_all`（証明済み）の対偶を取ると、
`l` が局所高さ `v_K(q_E)` と互いに素なら `l ∣ k` であり、
したがって安定な直線は `𝔽_l(1) = μ_l` に対応する。

☆機構は 1 行である: `k = l·m` なら `(x·Q^{-m})^l = 1`。 -/
theorem exists_rootOfUnity_mk_eq (Q : Kˣ) {l : ℕ} (hl : 0 < l)
    (x : Kˣ) (k : ℤ) (hxk : x ^ l = Q ^ k) (hdvd : (l : ℤ) ∣ k) :
    ∃ ζ : Kˣ, ζ ^ l = 1 ∧
      (QuotientGroup.mk x : Kˣ ⧸ Subgroup.zpowers Q) = QuotientGroup.mk ζ := by
  obtain ⟨m, rfl⟩ := hdvd
  refine ⟨x * (Q ^ m)⁻¹, ?_, ?_⟩
  · -- ★`(x·Q^{-m})^l = x^l·Q^{-lm} = 1`
    have hQ : ((Q ^ m)⁻¹ : Kˣ) ^ l = (Q ^ ((l : ℤ) * m))⁻¹ := by
      rw [inv_pow, ← zpow_natCast (Q ^ m) l, ← zpow_mul, mul_comm]
    rw [mul_pow, hxk, hQ, mul_inv_cancel]
  · -- ★両者の比は `Q^m` なので同じ類
    refine (QuotientGroup.eq (s := Subgroup.zpowers Q)).2 ?_
    have hval : x⁻¹ * (x * (Q ^ m)⁻¹) = (Q ^ m)⁻¹ := by
      rw [← mul_assoc, inv_mul_cancel, one_mul]
    rw [hval]
    exact Subgroup.inv_mem _ (Subgroup.zpow_mem _ (Subgroup.mem_zpowers Q) m)

/-! ## ★出典の紐付け(`.src`) -/

def exists_rootOfUnity_mk_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(l ∣ k なら類は μ_l の元の類である。★無条件)",
    sectionId := "genell-lemma-3-2" }

def mk_pow_injOn.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—μ_l の元は Kˣ/Q^ℤ の中でも区別される。★無条件)",
    sectionId := "genell-def-3-3" }

def normRep_of_pow_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—1 の冪根の類の代表元はそれ自身。★無条件)",
    sectionId := "genell-def-3-3" }

def tatePhi_of_pow_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—μ_l の点は tatePtPair ζⁱ (qζ⁻ⁱ)。★無条件)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
