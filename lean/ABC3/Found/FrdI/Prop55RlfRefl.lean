/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24SuppElt
import ABC3.Found.FrdI.Def24Rlf
import ABC3.Found.FrdI.Prop55RlfNd

/-!
# [FrdI] Proposition 5.5, (iii) —— `toSc` は `⪯` を(primary な相手について)反射する

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★`Prop55RlfNd.lean` は `hnd`(係数拡大した `Φ` が non-dilating)を
単系の 2 条 —— `hprim`(準素元の保存)と **`hrefl`(`⪯` の反射)** —— へ還元した。
本ファイルは **`hrefl` を perf-factorial の仮定から閉じる**。

## ★★★測って分かったこと(2026-08-25)—— `rlf-agree` は要らなかった

`Prop55RlfNd.lean` の ★3 は「perf-factorial の座標は**因子模型・錐模型**の上にあるので、
テンソル模型 `ScT ℝ≥0 Φ` との**一致**(節点 `rlf-agree`)を経由しないと使えない」と
書いていた。★**これは強すぎる見立てだった。**

要るのは**一致(同型)ではなく片道の準同型 1 本**である:

```
scToFactor : ℝ≥0 ⊗_ℕ M →+ (Prime M → ℝ≥0)      (r ⊗ m ↦ r • factorMap ι (Pf.of m))
```

これは `TensorProduct.liftAddHom` で作れる(`Def24PfT.lean` の `pfTToPf` と同じ型)。
★`MPrec (toSc x) (toSc y)` に `scToFactor` を当てれば
`μ x + (非負) = n • μ y` になり、`ℝ≥0` は
「和が `0` なら各項が `0`」なので **台の包含 `SuppElt x ⊆ SuppElt y`** が出る。
★★あとは在庫の `mprec_of_suppElt_eq_singleton`(同じ 1 点の台を持つ 2 元は `⪯`-同値)で
`MPrec x y` に戻る —— **相手が primary ならこれで足りる**。

★`Prop55RlfNd.lean` の `hrefl` は **primary な相手にしか使っていない**ので、
本ファイルの形で十分である。

## ★★★★残る `hprim` の測定(2026-08-25)—— こちらは**両向き**が要る

`hrefl` が片道で足りたのに対し、残る `hprim`
(`b` が準素元なら `toSc b` に `≼`-同値な準素元がある)は**そうはいかない**。

`IsPrimaryElt a = a ≠ 0 ∧ ∀ z ≠ 0, z ⪯ a → a ⪯ z` なので、後半には

* `z ≠ 0 → scToFactor z ≠ 0`(★**片道の橋の単射性**)
* 「`p` 成分が正なら `toSc b ⪯ z`」を示すための **`c` の構成**(★**逆向き**)

の**両方**が要る。★したがって `hprim` は本当に節点 `rlf-agree`
(`ℝ≥0 ⊗_ℕ M` と素点ごとの座標の**一致**)の側にある。

★★**弱めておいた**: `Prop55RlfNd.lean` の `hprim` は
「`toSc b` **そのもの**が準素元」ではなく
「`toSc b` に **`≼`-同値な準素元がある**」で足りる
(`scMap f (toSc b) ⪯ scMap f z ⪯ z ⪯ toSc b` と繋がるため)。
強い形からは `hprimWeak_of_isPrimaryElt` で降りる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal TensorProduct

universe v u w

section Refl

variable {M : Type w} [AddCommMonoid M] {ι : Prime M → Pf M → ℝ≥0}

/-! ## ★1. 台の写像 `μ : M →+ (Prime M → ℝ≥0)` -/

/-- ★★**元の因子分解** `μ a = factorMap ι (Pf.of a)`。 -/
noncomputable def suppFactorHom (H : IsPerfFactorialWith M ι) : M →+ RlfFactor M where
  toFun a := factorMap ι (Pf.of a)
  map_zero' := by
    show factorMap ι (Pf.of (0 : M)) = 0
    rw [map_zero]
    exact factorMap_zero H
  map_add' a b := by
    show factorMap ι (Pf.of (a + b)) = factorMap ι (Pf.of a) + factorMap ι (Pf.of b)
    rw [map_add]
    exact H.factorAdd _ _

@[simp] theorem suppFactorHom_apply (H : IsPerfFactorialWith M ι) (a : M) :
    suppFactorHom H a = factorMap ι (Pf.of a) := rfl

/-! ## ★2. `ℝ≥0 ⊗_ℕ M → (Prime M → ℝ≥0)` —— **片道の橋** -/

/-- ★`r ↦ (a ↦ r • μ a)`。 -/
noncomputable def rsmulFactorHom (H : IsPerfFactorialWith M ι) :
    ℝ≥0 →+ (M →+ RlfFactor M) where
  toFun r := r • suppFactorHom H
  map_zero' := by ext a; simp
  map_add' r r' := by ext a; simp [add_smul]

/-- ★★★★★**片道の橋** `ℝ≥0 ⊗_ℕ M →+ (Prime M → ℝ≥0)`。

★★**`Φ^rlf` の 2 模型の一致(`rlf-agree`)は要らない** —— 片道で足りる。 -/
noncomputable def scToFactor (H : IsPerfFactorialWith M ι) :
    ScT ℝ≥0 M →+ RlfFactor M :=
  TensorProduct.liftAddHom (R := ℕ) (rsmulFactorHom H) (by
    intro n r a
    show (n • r) • (suppFactorHom H a) = r • (suppFactorHom H (n • a))
    rw [map_nsmul, smul_comm r n, smul_assoc])

@[simp] theorem scToFactor_tmul (H : IsPerfFactorialWith M ι) (r : ℝ≥0) (a : M) :
    scToFactor H (r ⊗ₜ[ℕ] a) = r • factorMap ι (Pf.of a) := rfl

theorem scToFactor_toSc (H : IsPerfFactorialWith M ι) (a : M) :
    scToFactor H (toSc (S := ℝ≥0) a) = factorMap ι (Pf.of a) := by
  show scToFactor H ((1 : ℝ≥0) ⊗ₜ[ℕ] a) = _
  rw [scToFactor_tmul, one_smul]

/-! ## ★3. `⪯` から台の包含へ -/

/-- ★★★★**`ℝ≥0 ⊗` の上の `⪯` は台の包含を与える**。

★`ℝ≥0` は「和が `0` なら各項が `0`」なので、`μ x + (非負) = n • μ y` から
`μ y p = 0 → μ x p = 0` が出る。 -/
theorem suppElt_subset_of_mprec_sc (H : IsPerfFactorialWith M ι) {x y : M}
    (h : MPrec (toSc (S := ℝ≥0) x) (toSc y)) : SuppElt ι x ⊆ SuppElt ι y := by
  obtain ⟨n, -, c, hc⟩ := h
  have hkey := congrArg (scToFactor H) hc
  rw [map_add, map_nsmul, scToFactor_toSc, scToFactor_toSc] at hkey
  intro p hp
  by_contra hq
  have hy0 : factorMap ι (Pf.of y) p = 0 := by
    simpa [SuppElt, Supp] using hq
  have hev := congrFun hkey p
  rw [Pi.add_apply, Pi.smul_apply, hy0, smul_zero] at hev
  have hx0 : factorMap ι (Pf.of x) p = 0 :=
    (add_eq_zero.mp hev).1
  exact hp (by simpa [SuppElt, Supp] using hx0)

/-! ## ★4. `hrefl`(相手が primary の場合) -/

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii) の残り 1 条** ——
**`toSc` は `⪯` を(相手が準素元のとき)反射する**。

★これで `Prop55RlfNd.lean` の `hrefl` が閉じる。
★★`rlf-agree`(`Φ^rlf` の 2 模型の一致)は**要らなかった** —— 片道の橋で足りる。 -/
theorem mprec_of_mprec_sc_primary (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {b : M} {Pb : Prime M} (hbs : SuppElt ι b = {Pb}) (x : M)
    (h : MPrec (toSc (S := ℝ≥0) x) (toSc b)) : MPrec x b := by
  by_cases hx : x = 0
  · subst hx
    exact ⟨1, one_pos, b, by rw [one_smul, zero_add]⟩
  · have hsub := suppElt_subset_of_mprec_sc H h
    have hne : SuppElt ι x ≠ ∅ := suppElt_ne_empty H hdiv hx
    have hxs : SuppElt ι x = {Pb} := by
      rcases Set.subset_singleton_iff_eq.mp (hbs ▸ hsub) with h0 | h1
      · exact absurd h0 hne
      · exact h1
    exact mprec_of_suppElt_eq_singleton H hdiv hxs hbs

/-- ★★**`hperf` を使う版**(在庫の「準素元の台は 1 点」を経由する形)。

★★`hperf`(`M` が perfect)は**原文の常備仮定ではない** ——
台が 1 点であることを言う在庫の補題が要求するだけである。
上の `mprec_of_mprec_sc_primary` は**台が 1 点であること**しか使わないので、
そちらを直に与えられるならこの版は要らない。 -/
theorem mprec_of_mprec_sc_primary_of_perfect (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M)
    {b : M} (hb : IsPrimaryElt b) (x : M)
    (h : MPrec (toSc (S := ℝ≥0) x) (toSc b)) : MPrec x b :=
  mprec_of_mprec_sc_primary H hdiv (suppElt_eq_singleton_toPrime H hperf hdiv hb) x h

end Refl

/-! ## ★5. `MonoidOn` の層 —— `hrefl` を消した形 -/

section ReflOn

variable {D : Type u} [Category.{v} D]

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** ——
**`Φ` が perf-factorial なら、係数を `ℝ≥0` に拡げても non-dilating 性は保たれる**
(`hrefl` を消した形)。

★`Prop55RlfNd.lean` の 2 条のうち **`hrefl` はこれで消えた**。
残るのは `hprim`(準素元が `toSc` で準素元へ移ること)だけである。

★★`hperf`(`Φ(A)` が perfect)は、台が 1 点であることを言う在庫の補題
(`suppElt_eq_singleton_toPrime` / `suppElt_singleton_of_primary`)が要求するもので、
**原文の常備仮定ではない**。
★核心の `mprec_of_mprec_sc_primary` は **`hperf` を使わず「台が 1 点」だけを使う**形に
してあるので、台が 1 点であることを別の道で与えられればこの仮定は消える。 -/
theorem MonoidOn.scOn_isNonDilatingOn_of_perfFactorial (Φ : MonoidOn.{v, u, w} D)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := ℝ≥0) (Φ.map α)))
    (hsS : ∀ A : D, IsSharp (ScT ℝ≥0 (Φ.val A)))
    (hpf : ∀ A : D, IsPerfFactorial (Φ.val A))
    (hperf : ∀ A : D, IsPerfectMonoid (Φ.val A))
    (hprim : ∀ (A : D) (b : Φ.val A), IsPrimaryElt b →
      ∃ z : ScT ℝ≥0 (Φ.val A), IsPrimaryElt z
        ∧ MPrec (toSc (S := ℝ≥0) b) z ∧ MPrec z (toSc (S := ℝ≥0) b))
    (h : Φ.IsNonDilatingOn) : (phiScOn ℝ≥0 Φ hcharInj).IsNonDilatingOn := by
  refine MonoidOn.scOn_isNonDilatingOn Φ hcharInj
    (fun A => ((hpf A).choose_spec.divisorial).2) hsS hprim ?_ h
  intro A b hb x hx
  exact mprec_of_mprec_sc_primary_of_perfect (hpf A).choose_spec (hperf A)
    (hpf A).choose_spec.divisorial hb x hx

end ReflOn

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (iii)` の「`toSc` が `⪯` を反射する」。 -/
def mprec_of_mprec_sc_primary.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — toSc は ⪯ を反射する(相手が準素元のとき)",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★locator —— テンソル模型から因子模型への**片道の橋**
(`Definition 2.4` の座標をテンソル模型で使うのに要る。★2 模型の一致は不要)。 -/
def scToFactor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 48,
    item := "Definition 2.4 — ℝ≥0 ⊗_ℕ M から素点ごとの座標への準同型",
    sectionId := "frdi-def-2-4" }

end ABC3.Found.FrdI
