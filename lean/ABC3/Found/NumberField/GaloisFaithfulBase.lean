/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.GaloisFaithful

/-!
# 素点への作用は忠実(★**底が一般の数体**の版)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.115):
> since any automorphism of a number field that fixes all of the valuations of the
> number field is clearly equal to the identity automorphism

## ★★なぜ底を一般化するのか

`GaloisFaithful.lean` は同じ議論を**底 `ℚ`** で書いており、`[IsGalois ℚ K]` を要求する。
★しかし `Theorem 6.4, (i)` の `𝒟 = B(G)⁰` の対象 `L` は
**`F` の有限拡大**であって、**ℚ 上 Galois とは限らない**。

★★そこで「全素点を固定する自己同型 `τ` は恒等」を示すには

    M := L^{⟨τ⟩}

と置いて **`L/M` を巡回 Galois 拡大として見る**のが道筋である
(節点 `thm64-i-divslim` / `nf-spl-base`)。★本ファイルはその
**底が一般の数体の場合**の 3 本を用意する。

## ★★★中身は `ℚ` の場合と同じ

| 段 | 根拠 |
|---|---|
| `g·(e·f) = #Gal` | `Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn` |
| 完全分解 ⟹ `e·f = 1` | 上の等式と `g = [L:M]` |
| `#(stabilizer) = e·f` | `Ideal.card_stabilizer_eq`(mathlib、`IsGaloisGroup G R S` だけ要求) |

★`IsGaloisGroup (L ≃ₐ[M] L) (𝓞 M) (𝓞 L)` は **mathlib の instance がそのまま効く**
(底が `ℚ` のときのように `𝓞 ℚ ≅ ℤ` で移す必要が無い)。

## ☆残るもの —— 節点 `nf-spl-base`

本ファイルが要求する `hsplit`(`q` が `L` で完全分解する)の**存在**は、
底が `ℚ` なら在庫の `infinite_splitsCompletely`(★Chebotarev を使わない)が与える。
★一般の底 `M` では

1. `L̃` を `L/ℚ` の Galois 閉包とし、`L̃` で完全分解する有理素数 `p` を取る、
2. `p` は `L` でも完全分解する(`e`・`f` の乗法性)、
3. `P | p`、`q := P ∩ 𝓞 M` とすると `e(P/q) = f(P/q) = 1`

という降下が要る。★そこが節点 `nf-spl-base` である。

## ★実装上の注意(在庫と同じ)

* `Ideal.Quotient.field` は mathlib で**大域 instance ではない**ので
  `attribute [local instance]` が要る。
* `Ideal` への群作用は **`Pointwise` スコープ**にある(`open scoped Pointwise`)。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal
open scoped _root_.NumberField Pointwise

-- ★`Ideal.Quotient.field` は mathlib で**大域 instance ではない**(local instance)。
attribute [local instance] Ideal.Quotient.field

variable {M L : Type} [Field M] [Field L] [NumberField M] [NumberField L]
  [Algebra M L] [IsGalois M L]

/-! ## ★1. 完全分解なら `e·f = 1` -/

/-- ★★**完全分解なら `e·f = 1`**(底が一般の数体の版)。

★Galois の基本等式 `g·e·f = #Gal(L/M) = [L:M]` を当てるだけ。 -/
theorem ramificationIdxIn_mul_inertiaDegIn_eq_one_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hq : q ≠ ⊥)
    (hsplit : (q.primesOver (𝓞 L)).ncard = Module.finrank M L) :
    q.ramificationIdxIn (𝓞 L) * q.inertiaDegIn (𝓞 L) = 1 := by
  have hG := Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn
    (p := q) hq (𝓞 L) (L ≃ₐ[M] L)
  rw [IsGalois.card_aut_eq_finrank M L, hsplit] at hG
  have hpos : 0 < Module.finrank M L := Module.finrank_pos
  refine Nat.eq_of_mul_eq_mul_left hpos ?_
  rw [mul_one]
  exact hG

/-! ## ★2. 剰余体拡大は分離的 -/

omit [IsGalois M L] in
/-- ★★**剰余体拡大は分離的**(底が一般の数体の版)—— 剰余体は有限、有限体は完全だから。 -/
theorem isSeparable_residue_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hqne : q ≠ ⊥)
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.LiesOver q] (hPne : P ≠ ⊥) :
    Algebra.IsSeparable (𝓞 M ⧸ q) (𝓞 L ⧸ P) := by
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hPne
  haveI : Finite (𝓞 M ⧸ q) := Ideal.finiteQuotientOfFreeOfNeBot q hqne
  haveI : Finite (𝓞 L ⧸ P) := Ideal.finiteQuotientOfFreeOfNeBot P hPne
  haveI : PerfectField (𝓞 M ⧸ q) := PerfectField.ofFinite
  haveI : Algebra.IsAlgebraic (𝓞 M ⧸ q) (𝓞 L ⧸ P) := Algebra.IsAlgebraic.of_finite _ _
  infer_instance

/-! ## ★3. 分解群は自明、したがって素点を固定する自己同型は恒等 -/

/-- ★★★**完全分解する素点の分解群は自明**(底が一般の数体の版)。 -/
theorem subsingleton_stabilizer_of_splitsCompletely_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hqne : q ≠ ⊥)
    (hsplit : (q.primesOver (𝓞 L)).ncard = Module.finrank M L)
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.LiesOver q] (hPne : P ≠ ⊥) :
    Nat.card (MulAction.stabilizer (L ≃ₐ[M] L) P) = 1 := by
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hPne
  haveI := isSeparable_residue_base q hqne P hPne
  rw [Ideal.card_stabilizer_eq q hqne P]
  exact ramificationIdxIn_mul_inertiaDegIn_eq_one_base q hqne hsplit

/-- ★★★★**完全分解する素点を 1 つ固定する自己同型は恒等**(底が一般の数体の版)。

★★これが `Theorem 6.4, (i)` の「数体の自己同型で全素点を固定するものは恒等」の
**最後の段**である —— `τ ∈ Aut(L/F)` が全素点を固定するとき
`M := L^{⟨τ⟩}` と置いて本定理を当てればよい。
★残るのは `hsplit`(`q` が `L` で完全分解する)の**存在**だけで、
それが節点 `nf-spl-base` である。 -/
theorem eq_one_of_fixes_prime_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hqne : q ≠ ⊥)
    (hsplit : (q.primesOver (𝓞 L)).ncard = Module.finrank M L)
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.LiesOver q] (hPne : P ≠ ⊥)
    (σ : L ≃ₐ[M] L) (hσ : σ • P = P) : σ = 1 := by
  have hmem : σ ∈ MulAction.stabilizer (L ≃ₐ[M] L) P := hσ
  have hcard := subsingleton_stabilizer_of_splitsCompletely_base q hqne hsplit P hPne
  haveI : Subsingleton (MulAction.stabilizer (L ≃ₐ[M] L) P) :=
    Nat.card_eq_one_iff_unique.mp hcard |>.1
  exact congrArg Subtype.val
    (Subsingleton.elim (⟨σ, hmem⟩ : MulAction.stabilizer (L ≃ₐ[M] L) P) 1)

/-! ## ★4. 完全分解は塔を降りる

★★節点 `nf-spl-base` の第 2 段。`M ⊆ L ⊆ N` で `N` の側の素点が `e = f = 1` なら、
**間の `L` でも `q` は完全分解する**。

★中身は mathlib の**塔での `e`・`f` の乗法性**
(`Ideal.ramificationIdx_algebra_tower'` / `Ideal.inertiaDeg_algebra_tower`)と、
Galois の基本等式 `g·e·f = #Gal(L/M) = [L:M]` だけである。 -/

section Tower

variable {M L N : Type} [Field M] [Field L] [Field N]
  [NumberField M] [NumberField L] [NumberField N]
  [Algebra M L] [Algebra L N] [Algebra M N] [IsScalarTower M L N]
  [IsGalois M L]

/-- ★★★★**完全分解は塔を降りる** —— `q` の上の `N` の素点 `R` が `e = f = 1` なら、
`q` は `L` で完全分解する。 -/
theorem splitsCompletely_of_tower
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hq : q ≠ ⊥)
    (R : Ideal (𝓞 N)) [R.IsPrime]
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.IsMaximal] [P.LiesOver q] [R.LiesOver P]
    (he : q.ramificationIdx R = 1) (hf : q.inertiaDeg R = 1) :
    (q.primesOver (𝓞 L)).ncard = Module.finrank M L := by
  haveI := RingOfIntegers.inst_isScalarTower M L N
  have heP : q.ramificationIdx P = 1 := by
    have htower := Ideal.ramificationIdx_algebra_tower' q P R
    rw [he] at htower
    exact Nat.eq_one_of_mul_eq_one_right htower.symm
  have hfP : q.inertiaDeg P = 1 := by
    have htower := Ideal.inertiaDeg_algebra_tower q P R
    rw [hf] at htower
    exact Nat.eq_one_of_mul_eq_one_right htower.symm
  have hG := Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn
    (p := q) hq (𝓞 L) (L ≃ₐ[M] L)
  rw [Ideal.ramificationIdxIn_eq_ramificationIdx q P (L ≃ₐ[M] L),
    Ideal.inertiaDegIn_eq_inertiaDeg q P (L ≃ₐ[M] L), heP, hfP,
    mul_one, mul_one, IsGalois.card_aut_eq_finrank M L] at hG
  exact hG

/-- ★★★★★**塔の上で完全分解する素点があれば、その下の素点を固定する自己同型は恒等**。

★★これが `Theorem 6.4, (i)` の「数体の自己同型で全素点を固定するものは恒等」に
そのまま当たる形である —— `τ ∈ Aut(L/F)` が全素点を固定するとき

    M := L^{⟨τ⟩}、  N := L/ℚ の Galois 閉包

と取れば、`N` で完全分解する有理素数(在庫の `infinite_splitsCompletely`、
★Chebotarev を使わない)から `he` / `hf` が出る。 -/
theorem eq_one_of_fixes_prime_tower
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hq : q ≠ ⊥)
    (R : Ideal (𝓞 N)) [R.IsPrime]
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.IsMaximal] [P.LiesOver q] [R.LiesOver P]
    (hPne : P ≠ ⊥)
    (he : q.ramificationIdx R = 1) (hf : q.inertiaDeg R = 1)
    (σ : L ≃ₐ[M] L) (hσ : σ • P = P) : σ = 1 :=
  eq_one_of_fixes_prime_base q hq
    (splitsCompletely_of_tower q hq R P he hf) P hPne σ hσ

end Tower

/-! ## ★5. ℚ の側の完全分解から降ろす

★★節点 `nf-spl-base` の第 3 段。`N/ℚ` が Galois で有理素数 `p` が `N` で完全分解するなら、
中間体 `M` の下の素点 `q := R ∩ 𝓞 M` でも `e = f = 1` である。

★中身は**塔での `e`・`f` の乗法性をもう一度**当てるだけである
(`ℤ ⊆ 𝓞 M ⊆ 𝓞 N`)。 -/

section FromQ

/-- ★★★**ℚ の側の完全分解から、中間体 `M` の下の素点での `e = f = 1` を出す**。 -/
theorem ramificationIdx_inertiaDeg_eq_one_of_splitsCompletely
    {M N : Type} [Field M] [Field N] [NumberField M] [NumberField N]
    [Algebra M N] [IsGalois ℚ N] {p : ℕ} [Fact p.Prime]
    (hsplit : (Ideal.primesOver (Ideal.span {(p : ℤ)}) (𝓞 N)).ncard = Module.finrank ℚ N)
    (R : Ideal (𝓞 N)) [R.IsPrime] [R.LiesOver (Ideal.span {(p : ℤ)})] (hRne : R ≠ ⊥) :
    (R.under (𝓞 M)).ramificationIdx R = 1 ∧ (R.under (𝓞 M)).inertiaDeg R = 1 := by
  haveI := RingOfIntegers.inst_isScalarTower ℚ M N
  haveI : IsScalarTower ℤ (𝓞 M) (𝓞 N) := inferInstance
  haveI : R.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hRne
  haveI : (Ideal.span {(p : ℤ)} : Ideal ℤ).IsMaximal := Int.ideal_span_isMaximal_of_prime p
  haveI : (R.under (𝓞 M)).IsMaximal := isMaximal_comap_of_isIntegral_of_isMaximal R
  have hef := ramificationIdxIn_mul_inertiaDegIn_eq_one (K := N) hsplit
  have heZ : (Ideal.span {(p : ℤ)} : Ideal ℤ).ramificationIdx R = 1 := by
    rw [← Ideal.ramificationIdxIn_eq_ramificationIdx _ R (N ≃ₐ[ℚ] N)]
    exact Nat.eq_one_of_mul_eq_one_right hef
  have hfZ : (Ideal.span {(p : ℤ)} : Ideal ℤ).inertiaDeg R = 1 := by
    rw [← Ideal.inertiaDegIn_eq_inertiaDeg _ R (N ≃ₐ[ℚ] N)]
    exact Nat.eq_one_of_mul_eq_one_left hef
  constructor
  · have ht := Ideal.ramificationIdx_algebra_tower' (Ideal.span {(p : ℤ)}) (R.under (𝓞 M)) R
    rw [heZ] at ht
    exact Nat.eq_one_of_mul_eq_one_left ht.symm
  · have ht := Ideal.inertiaDeg_algebra_tower (Ideal.span {(p : ℤ)}) (R.under (𝓞 M)) R
    rw [hfZ] at ht
    exact Nat.eq_one_of_mul_eq_one_left ht.symm

/-- ★★★★★★**数体の自己同型で素点を固定するものは恒等**(完全分解する有理素数を経由)。

★★★これが `Theorem 6.4, (i)` の「clearly」の中身である。使うときは

    M := L^{⟨τ⟩}(`τ` の固定体)、 N := `L/ℚ` の Galois 閉包

と取ればよい。★`hsplit`(`N` で完全分解する有理素数の存在)は在庫の
`infinite_splitsCompletely`(**Chebotarev を使わない**)が与える。 -/
theorem eq_one_of_fixes_prime_of_splitsCompletely
    {M L N : Type} [Field M] [Field L] [Field N]
    [NumberField M] [NumberField L] [NumberField N]
    [Algebra M L] [Algebra L N] [Algebra M N] [IsScalarTower M L N]
    [IsGalois M L] [IsGalois ℚ N] {p : ℕ} [Fact p.Prime]
    (hsplit : (Ideal.primesOver (Ideal.span {(p : ℤ)}) (𝓞 N)).ncard = Module.finrank ℚ N)
    (R : Ideal (𝓞 N)) [R.IsPrime] [R.LiesOver (Ideal.span {(p : ℤ)})] (hRne : R ≠ ⊥)
    (σ : L ≃ₐ[M] L) (hσ : σ • (R.under (𝓞 L)) = R.under (𝓞 L)) : σ = 1 := by
  haveI := RingOfIntegers.inst_isScalarTower M L N
  haveI := RingOfIntegers.inst_isScalarTower ℚ M N
  haveI hZL : IsScalarTower ℤ (𝓞 L) (𝓞 N) := inferInstance
  haveI hZM : IsScalarTower ℤ (𝓞 M) (𝓞 N) := inferInstance
  haveI : R.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hRne
  have hspan : (Ideal.span {(p : ℤ)} : Ideal ℤ) ≠ ⊥ := by simp [NeZero.ne p]
  set P : Ideal (𝓞 L) := R.under (𝓞 L) with hPdef
  haveI : P.IsMaximal := isMaximal_comap_of_isIntegral_of_isMaximal R
  haveI : R.LiesOver P := inferInstance
  haveI hPo : P.LiesOver (Ideal.span {(p : ℤ)}) :=
    ⟨by rw [hPdef, Ideal.under_under]; exact Ideal.LiesOver.over⟩
  have hPne : P ≠ ⊥ := Ideal.ne_bot_of_liesOver_of_ne_bot hspan P
  set q : Ideal (𝓞 M) := P.under (𝓞 M) with hqdef
  haveI : q.IsMaximal := isMaximal_comap_of_isIntegral_of_isMaximal P
  haveI : P.LiesOver q := inferInstance
  haveI hqo : q.LiesOver (Ideal.span {(p : ℤ)}) :=
    ⟨by rw [hqdef, hPdef, Ideal.under_under, Ideal.under_under]; exact Ideal.LiesOver.over⟩
  have hqne : q ≠ ⊥ := Ideal.ne_bot_of_liesOver_of_ne_bot hspan q
  have hqR : q = R.under (𝓞 M) := by rw [hqdef, hPdef, Ideal.under_under]
  obtain ⟨he, hf⟩ := ramificationIdx_inertiaDeg_eq_one_of_splitsCompletely (M := M) hsplit R hRne
  rw [← hqR] at he hf
  exact eq_one_of_fixes_prime_tower q hqne R P hPne he hf σ hσ

end FromQ

/-! ## ★6. 固定体は Galois(Artin)

★★節点 `nf-spl-base` の使い方の第 1 段 —— `τ` の固定体 `M := L^{⟨τ⟩}` を取ると
`L/M` は Galois で `Gal(L/M) ≃ ⟨τ⟩` になる。 -/

open IntermediateField in
/-- ★★★**部分群の固定体の上で `L` は Galois**(Artin の定理、数体の場合)。

★mathlib は `IntermediateField.finrank_fixedField_eq_card`(`[L : L^G] = #G`)と
`IntermediateField.subgroupEquivAlgEquiv`(`G ≃* Gal(L/L^G)`)を持つので、
`IsGalois.of_card_aut_eq_finrank` に流すだけである。 -/
theorem isGalois_fixedField (L : Type) [Field L] [NumberField L]
    (G : Subgroup (L ≃ₐ[ℚ] L)) :
    IsGalois (IntermediateField.fixedField G) L := by
  refine IsGalois.of_card_aut_eq_finrank (↥(fixedField G)) L ?_
  rw [IntermediateField.finrank_fixedField_eq_card G]
  exact Nat.card_congr (IntermediateField.subgroupEquivAlgEquiv G).symm.toEquiv

/-! ## ★7. 配線 —— `τ` から始めて `τ = 1` まで

★★★これが `Theorem 6.4, (i)` の「clearly」の**完成形**である。
`τ` の固定体を取り(★`isGalois_fixedField`)、
`L` を含む ℚ 上 Galois な `N` で完全分解する有理素数を経由して(★§5)、
最後に `subgroupEquivAlgEquiv` の単射性で `τ` に戻す。 -/

open IntermediateField in
/-- ★★★★★★★**数体の自己同型で素点を固定するものは恒等**(配線まで込み)。

★入力は「`L` を含む ℚ 上 Galois な `N` と、そこで完全分解する有理素数」だけである。
★★`N` は `L/ℚ` の Galois 閉包を取ればよく、完全分解する有理素数の存在は
在庫の `infinite_splitsCompletely`(**Chebotarev を使わない**)が与える。

★★★中身:

| 段 | 根拠 |
|---|---|
| `M := L^{⟨τ⟩}` の上で `L` は Galois | `isGalois_fixedField`(Artin) |
| `σ := τ` を `Gal(L/M)` の元として見る | `IntermediateField.subgroupEquivAlgEquiv` |
| `σ` は `R ∩ 𝓞 L` を固定する | `τ` と同じ底写像なので `rfl` |
| `σ = 1` | `eq_one_of_fixes_prime_of_splitsCompletely` |
| `τ = 1` | `subgroupEquivAlgEquiv` の単射性 | -/
theorem eq_one_of_fixes_prime_of_galoisClosure
    {L N : Type} [Field L] [Field N] [NumberField L] [NumberField N]
    [Algebra L N] [IsGalois ℚ N]
    {p : ℕ} [Fact p.Prime]
    (hsplit : (Ideal.primesOver (Ideal.span {(p : ℤ)}) (𝓞 N)).ncard = Module.finrank ℚ N)
    (R : Ideal (𝓞 N)) [R.IsPrime] [R.LiesOver (Ideal.span {(p : ℤ)})] (hRne : R ≠ ⊥)
    (τ : L ≃ₐ[ℚ] L) (hτ : τ • (R.under (𝓞 L)) = R.under (𝓞 L)) : τ = 1 := by
  haveI : IsScalarTower ℚ L N := IsScalarTower.of_algebraMap_eq (fun x => by
    show algebraMap ℚ N x = _
    exact (IsScalarTower.algebraMap_apply ℚ L N x))
  set G : Subgroup (L ≃ₐ[ℚ] L) := Subgroup.zpowers τ with hG
  haveI : IsGalois (IntermediateField.fixedField G) L := isGalois_fixedField L G
  set σ : L ≃ₐ[↥(IntermediateField.fixedField G)] L :=
    IntermediateField.subgroupEquivAlgEquiv G ⟨τ, Subgroup.mem_zpowers τ⟩ with hσdef
  have hσact : σ • (R.under (𝓞 L)) = R.under (𝓞 L) := hτ
  have hσ1 : σ = 1 :=
    eq_one_of_fixes_prime_of_splitsCompletely hsplit R hRne σ hσact
  have hsub := (IntermediateField.subgroupEquivAlgEquiv G).injective
    ((hσdef ▸ hσ1).trans (map_one _).symm)
  exact congrArg Subtype.val hsub

/-! ### ★出典の紐付け -/

def eq_one_of_fixes_prime_of_galoisClosure.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 数体の自己同型で素点を固定するものは恒等(配線まで)",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_prime_of_galoisClosure.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isGalois_fixedField(Artin。τ の固定体の上で L は Galois)"
      (.inProject "ABC3" "ABC3.Found.NF.isGalois_fixedField") 114,
    .citation "[ABC3]" "eq_one_of_fixes_prime_of_splitsCompletely(完全分解を経由する本体)"
      (.inProject "ABC3" "ABC3.Found.NF.eq_one_of_fixes_prime_of_splitsCompletely") 114,
    .citation "[mathlib]" "IntermediateField.subgroupEquivAlgEquiv(G ≃* Gal(L/L^G))"
      (.inMathlib "IntermediateField.subgroupEquivAlgEquiv") 114,
    .implicitStep
      ("★入力の N は L/ℚ の Galois 閉包を取ればよく、完全分解する有理素数の存在は" ++
       "在庫の infinite_splitsCompletely(Chebotarev を使わない)が与える") 116 ]

def isGalois_fixedField.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 固定体 L^{⟨τ⟩} の上で L は Galois(Artin)",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_prime_of_splitsCompletely.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 数体の自己同型で素点を固定するものは恒等",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_prime_of_splitsCompletely.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "infinite_splitsCompletely(完全分解する有理素数の存在。Chebotarev を使わない)"
      (.inProject "ABC3" "ABC3.Found.NF.infinite_splitsCompletely") 116,
    .citation "[mathlib]" "Ideal.ramificationIdx_algebra_tower' / Ideal.inertiaDeg_algebra_tower"
      (.inMathlib "Ideal.ramificationIdx_algebra_tower'") 116,
    .implicitStep
      ("★使うときは M := L^{⟨τ⟩}(τ の固定体)、N := L/ℚ の Galois 閉包 と取る。" ++
       "残るのはその 2 つを組む配線だけである") 114 ]

def splitsCompletely_of_tower.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (i) — 完全分解は塔を降りる",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_prime_tower.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 塔の上の完全分解から素点を固定する自己同型は恒等",
    sectionId := "frdi-thm-6-4" }

def splitsCompletely_of_tower.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.ramificationIdx_algebra_tower'(塔での e の乗法性)"
      (.inMathlib "Ideal.ramificationIdx_algebra_tower'") 116,
    .citation "[mathlib]" "Ideal.inertiaDeg_algebra_tower(塔での f の乗法性)"
      (.inMathlib "Ideal.inertiaDeg_algebra_tower") 116,
    .citation "[mathlib]" "Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn(基本等式)"
      (.inMathlib "Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn") 116 ]

/-- ★★★★locator —— `Theorem 6.4, (i)` の「全素点を固定する自己同型は恒等」
(★**底が一般の数体**の版。`ℚ` 上 Galois とは限らない `L` に当てるため)。 -/
def eq_one_of_fixes_prime_base.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 完全分解する素点を固定する自己同型は恒等(一般の底)",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_prime_base.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.card_stabilizer_eq(#stabilizer = e·f)"
      (.inMathlib "Ideal.card_stabilizer_eq") 114,
    .citation "[mathlib]" "IsGaloisGroup (L ≃ₐ[M] L) (𝓞 M) (𝓞 L)(instance)"
      (.inMathlib "IsGaloisGroup") 114,
    .citation "[ABC3]" "ℚ を底とする版(同じ議論)"
      (.inProject "ABC3" "ABC3.Found.NF.eq_one_of_fixes_prime") 114,
    .implicitStep
      ("★仮定 hsplit(q が L で完全分解する)の存在が残る。底が ℚ なら在庫の " ++
       "infinite_splitsCompletely(Chebotarev を使わない)が与えるが、" ++
       "一般の底 M では Galois 閉包経由の降下が要る ＝ 節点 nf-spl-base") 116 ]

end ABC3.Found.NF
