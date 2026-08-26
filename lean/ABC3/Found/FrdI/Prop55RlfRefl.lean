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

## ★★★★★★★残る `hprim` も閉じた(2026-08-25、第 2 弾)

一時「`hprim` は `rlf-agree`(2 模型の**一致**)の側にある」と書いたが、
★**これも強すぎる見立てだった**。要ったのは 3 つで、いずれも片道の橋で足りる:

1. **核が自明** —— `scToFactor z = 0 → z = 0`(`scToFactor_eq_zero`)。
   `ℝ≥0` と `Prime M → ℝ≥0` は**負の元を持たない**ので、
   `Σ rᵢ • μ(aᵢ) = 0` から各項が `0` になり、`rᵢ = 0` または `aᵢ = 0`。
2. **`0` でないテンソルには `0` でない単項式が下にある**(`exists_tmul_mle`)——
   テンソルの帰納法だけ。
3. **Archimedes 性** —— `r ⊗ a ≤ w` かつ `r ≠ 0` なら
   `N·r ≥ 1` となる `N` を取って `N • w = toSc a + (残り)`(`mprec_toSc_of_tmul_mle`)。

★★これで `IsPrimaryElt (toSc b)` が出る(`isPrimaryElt_toSc`)——
`w ≠ 0` から単項式 `r ⊗ a ≤ w` を取り、台が `b` と同じ 1 点になることを見て
在庫の `mprec_of_suppElt_eq_singleton` で `M` の中で `b ⪯ a`、
それを `toSc` で送り、Archimedes 性で `toSc a ⪯ w` に繋ぐ。

★★★あわせて `hsS`(`ScT ℝ≥0 M` が sharp)も核の自明性から出た(`isSharp_scT`)。

★★★★**結論**: `Proposition 5.5, (iii)` の `hnd` は
`Φ` が perf-factorial（＋各 `Φ(A)` が perfect）だけから出る
(`MonoidOn.scOn_isNonDilatingOn_of_perfFactorial`)。
節点 `rlf-agree`(2 模型の**一致**)は**どちらにも要らなかった**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal TensorProduct

universe v u w u2 v2

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

/-! ## ★2-b. 片道の橋は**核が自明** -/

/-- ★★`μ = factorMap ∘ Pf.of` は単射(`M` が divisorial かつ perf-factorial なら)。 -/
theorem suppFactorHom_injective (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M) :
    Function.Injective (suppFactorHom H) := by
  intro a b hab
  exact Pf.of_injective_of_divisorial hdiv (H.factorInj hab)

/-- ★★★★★**片道の橋の核は自明** —— `scToFactor z = 0` なら `z = 0`。

★★`ℝ≥0` と `Prime M → ℝ≥0` は**負の元を持たない**ので、
`Σ rᵢ • μ(aᵢ) = 0` から各項が `0` になり、`rᵢ = 0` または `aᵢ = 0`、
どちらでも `rᵢ ⊗ aᵢ = 0` になる。

★★★これが節点 `rlf-agree` の**片割れ**である ——
`Proposition 5.5, (iii)` の残り 1 条 `hprim` が要求する 2 つのうち、
「`z ≠ 0 → scToFactor z ≠ 0`」の側がこれで閉じる。 -/
theorem scToFactor_eq_zero (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {z : ScT ℝ≥0 M} (h : scToFactor H z = 0) : z = 0 := by
  induction z using TensorProduct.induction_on with
  | zero => rfl
  | tmul r a =>
    rw [scToFactor_tmul] at h
    by_cases hr : r = 0
    · rw [hr, TensorProduct.zero_tmul]
    · have ha : factorMap ι (Pf.of a) = 0 := by
        funext p
        have hp := congrFun h p
        rw [Pi.smul_apply, Pi.zero_apply, smul_eq_mul] at hp
        exact (mul_eq_zero.mp hp).resolve_left hr
      have : a = 0 := by
        refine suppFactorHom_injective H hdiv ?_
        rw [suppFactorHom_apply, suppFactorHom_apply, ha, map_zero]
        exact (factorMap_zero H).symm
      rw [this, TensorProduct.tmul_zero]
  | add x y hx hy =>
    rw [map_add] at h
    have h0 : ∀ p : Prime M, scToFactor H x p = 0 ∧ scToFactor H y p = 0 := by
      intro p
      have hp := congrFun h p
      rw [Pi.add_apply, Pi.zero_apply] at hp
      exact add_eq_zero.mp hp
    rw [hx (funext fun p => (h0 p).1), hy (funext fun p => (h0 p).2), add_zero]

/-- ★★★★同上、対偶の形。 -/
theorem scToFactor_ne_zero (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {z : ScT ℝ≥0 M} (h : z ≠ 0) : scToFactor H z ≠ 0 :=
  fun h0 => h (scToFactor_eq_zero H hdiv h0)

/-- ★★★★**`ScT ℝ≥0 M` は sharp**(`M` が perf-factorial かつ divisorial なら)。

★可逆元 `a` は `a + b = 0` を満たすので、`ℝ≥0` 側で `scToFactor a = 0`、
核が自明だから `a = 0`。 -/
theorem isSharp_scT (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M) :
    IsSharp (ScT ℝ≥0 M) := by
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h0 : ((u : ScT ℝ≥0 M) + (u.neg : ScT ℝ≥0 M)) = 0 := u.val_neg
  have hψ : scToFactor H (u : ScT ℝ≥0 M) + scToFactor H (u.neg : ScT ℝ≥0 M) = 0 := by
    rw [← map_add, h0, map_zero]
  refine scToFactor_eq_zero H hdiv (funext fun p => ?_)
  have hp := congrFun hψ p
  rw [Pi.add_apply, Pi.zero_apply] at hp
  exact (add_eq_zero.mp hp).1

/-- ★★★★**`toSc` は単射**(`M` が perf-factorial かつ divisorial なら)。

★`scToFactor ∘ toSc = μ` で `μ` が単射だから。 -/
theorem toSc_injective (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M) :
    Function.Injective (toSc (S := ℝ≥0) (M := M)) := by
  intro a b hab
  refine suppFactorHom_injective H hdiv ?_
  rw [suppFactorHom_apply, suppFactorHom_apply, ← scToFactor_toSc H, ← scToFactor_toSc H, hab]

/-- ★★★**`toSc` は `0` でない元を `0` でない元へ送る**(`M` が perf-factorial なら)。

★`hprim` の第 1 条(`toSc b ≠ 0`)がこれで閉じる。 -/
theorem toSc_ne_zero (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {b : M} (hb : b ≠ 0) : toSc (S := ℝ≥0) b ≠ 0 := by
  intro h0
  refine hb (suppFactorHom_injective H hdiv ?_)
  rw [suppFactorHom_apply, suppFactorHom_apply, ← scToFactor_toSc H, ← scToFactor_toSc H,
    map_zero, h0, map_zero]

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

/-! ## ★3-b. `hprim` —— `toSc b` は準素元 -/

/-- ★★★**`0` でないテンソルには `0` でない単項式が下にある**。 -/
theorem exists_tmul_mle {w : ScT ℝ≥0 M} (hw : w ≠ 0) :
    ∃ (r : ℝ≥0) (a : M), r ≠ 0 ∧ a ≠ 0 ∧ MLe (r ⊗ₜ[ℕ] a) w := by
  induction w using TensorProduct.induction_on with
  | zero => exact absurd rfl hw
  | tmul r a =>
    by_cases hr : r = 0
    · exact absurd (by rw [hr, TensorProduct.zero_tmul]) hw
    by_cases ha : a = 0
    · exact absurd (by rw [ha, TensorProduct.tmul_zero]) hw
    exact ⟨r, a, hr, ha, 0, add_zero _⟩
  | add x y hx hy =>
    by_cases hx0 : x = 0
    · obtain ⟨r, a, hr, ha, c, hc⟩ := hy (by rw [hx0, zero_add] at hw; exact hw)
      exact ⟨r, a, hr, ha, c, by rw [hx0, zero_add]; exact hc⟩
    · obtain ⟨r, a, hr, ha, c, hc⟩ := hx hx0
      exact ⟨r, a, hr, ha, c + y, by rw [← add_assoc, hc]⟩

/-- ★★★★**単項式が下にあれば `toSc a` は `⪯` で下にある**。

★`ℝ≥0` は Archimedes 的なので、`N·r ≥ 1` となる `N` を取れば
`N • w = toSc a + (残り)` と書ける。 -/
theorem mprec_toSc_of_tmul_mle {r : ℝ≥0} (hr : r ≠ 0) {a : M} {w : ScT ℝ≥0 M}
    (h : MLe (r ⊗ₜ[ℕ] a) w) : MPrec (toSc (S := ℝ≥0) a) w := by
  obtain ⟨d, hd⟩ := h
  obtain ⟨N, hN⟩ := exists_nat_ge (1 / r : ℝ≥0)
  have hrpos : (0 : ℝ≥0) < r := pos_iff_ne_zero.mpr hr
  set s : ℝ≥0 := ((N : ℝ≥0) + 1) * r with hs
  have hs1 : (1 : ℝ≥0) ≤ s := by
    have h1 : (1 / r : ℝ≥0) ≤ (N : ℝ≥0) + 1 := le_trans hN (le_add_of_nonneg_right zero_le')
    have := mul_le_mul_right' h1 r
    rwa [div_mul_cancel₀ _ hr] at this
  refine ⟨N + 1, Nat.succ_pos N, (s - 1) ⊗ₜ[ℕ] a + (N + 1) • d, ?_⟩
  have hsm : (N + 1) • (r ⊗ₜ[ℕ] a) = s ⊗ₜ[ℕ] a := by
    have hst : ((N + 1 : ℕ) • r) ⊗ₜ[ℕ] a = (N + 1 : ℕ) • (r ⊗ₜ[ℕ] a) :=
      TensorProduct.smul_tmul' (N + 1 : ℕ) r a
    rw [← hst, hs]
    congr 1
    rw [nsmul_eq_mul]
    push_cast
    ring
  have hsplit : s ⊗ₜ[ℕ] a = toSc (S := ℝ≥0) a + (s - 1) ⊗ₜ[ℕ] a := by
    show s ⊗ₜ[ℕ] a = (1 : ℝ≥0) ⊗ₜ[ℕ] a + (s - 1) ⊗ₜ[ℕ] a
    rw [← TensorProduct.add_tmul, add_tsub_cancel_of_le hs1]
  calc toSc (S := ℝ≥0) a + ((s - 1) ⊗ₜ[ℕ] a + (N + 1) • d)
      = (toSc (S := ℝ≥0) a + (s - 1) ⊗ₜ[ℕ] a) + (N + 1) • d := by rw [add_assoc]
    _ = (N + 1) • (r ⊗ₜ[ℕ] a) + (N + 1) • d := by rw [hsm, hsplit]
    _ = (N + 1) • (r ⊗ₜ[ℕ] a + d) := (smul_add _ _ _).symm
    _ = (N + 1) • w := by rw [hd]

/-- ★★★**単項式が `toSc b` の下にあれば台も含まれる**。 -/
theorem suppElt_subset_of_mprec_tmul (H : IsPerfFactorialWith M ι) {r : ℝ≥0} (hr : r ≠ 0)
    {a b : M} (h : MPrec (r ⊗ₜ[ℕ] a) (toSc (S := ℝ≥0) b)) : SuppElt ι a ⊆ SuppElt ι b := by
  obtain ⟨n, -, c, hc⟩ := h
  have hkey := congrArg (scToFactor H) hc
  rw [map_add, map_nsmul, scToFactor_tmul, scToFactor_toSc] at hkey
  intro p hp
  by_contra hq
  have hy0 : factorMap ι (Pf.of b) p = 0 := by simpa [SuppElt, Supp] using hq
  have hev := congrFun hkey p
  rw [Pi.add_apply, Pi.smul_apply, Pi.smul_apply, hy0, smul_zero] at hev
  have hx0 : r • factorMap ι (Pf.of a) p = 0 := (add_eq_zero.mp hev).1
  rw [smul_eq_mul] at hx0
  exact hp (by simpa [SuppElt, Supp] using (mul_eq_zero.mp hx0).resolve_left hr)

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii) の残り 1 条 `hprim`** ——
**`b` の台が 1 点なら `toSc b` は準素元**。

★★これで `Prop55RlfNd.lean` の 2 条(`hprim` / `hrefl`)が**両方とも閉じた**。

★道筋:
1. `toSc b ≠ 0` …… `scToFactor` の核が自明(`toSc_ne_zero`)
2. `w ≠ 0` から `0` でない単項式 `r ⊗ a ≤ w` を取る(`exists_tmul_mle`)
3. `SuppElt a = SuppElt b = {p}` を出し、在庫の
   `mprec_of_suppElt_eq_singleton` で `M` の中で `b ⪯ a`
4. `toSc` は `⪯` を保つので `toSc b ⪯ toSc a`
5. Archimedes 性で `toSc a ⪯ w`(`mprec_toSc_of_tmul_mle`) -/
theorem isPrimaryElt_toSc (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {b : M} {Pb : Prime M} (hbs : SuppElt ι b = {Pb}) (hb0 : b ≠ 0) :
    IsPrimaryElt (toSc (S := ℝ≥0) b) := by
  refine ⟨toSc_ne_zero H hdiv hb0, fun w hw hwb => ?_⟩
  obtain ⟨r, a, hr, ha, hle⟩ := exists_tmul_mle hw
  have hsub : SuppElt ι a ⊆ SuppElt ι b :=
    suppElt_subset_of_mprec_tmul H hr (mprec_trans (mprec_of_mle hle) hwb)
  have hne : SuppElt ι a ≠ ∅ := suppElt_ne_empty H hdiv ha
  have has : SuppElt ι a = {Pb} := by
    rcases Set.subset_singleton_iff_eq.mp (hbs ▸ hsub) with h0 | h1
    · exact absurd h0 hne
    · exact h1
  have hba : MPrec b a := mprec_of_suppElt_eq_singleton H hdiv hbs has
  exact mprec_trans (mprec_map (toSc (S := ℝ≥0)) hba) (mprec_toSc_of_tmul_mle hr hle)

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

/-- ★★★★★★★**[FrdI] Proposition 5.5, (iii)** ——
**`Φ` が perf-factorial なら、係数を `ℝ≥0` に拡げても non-dilating 性は保たれる**。

★★★**`Prop55RlfNd.lean` の 2 条(`hprim` / `hrefl`)は両方とも消えた。**
`hsS`(`ScT ℝ≥0 Φ` が sharp)も片道の橋の核が自明なことから出るので消えた。

★★残る `hperf`(`Φ(A)` が perfect)は、**準素元の台が 1 点**であることを言う
在庫の補題(`suppElt_eq_singleton_toPrime` / `suppElt_singleton_of_primary`)が
要求するものであって、**原文の常備仮定ではない**。
★核心の `isPrimaryElt_toSc` / `mprec_of_mprec_sc_primary` は
**`hperf` を使わず「台が 1 点」だけを使う**形にしてあるので、
台が 1 点であることを別の道で与えられればこの仮定も消える。 -/
theorem MonoidOn.scOn_isNonDilatingOn_of_perfFactorial (Φ : MonoidOn.{v, u, w} D)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := ℝ≥0) (Φ.map α)))
    (hpf : ∀ A : D, IsPerfFactorial (Φ.val A))
    (hperf : ∀ A : D, IsPerfectMonoid (Φ.val A))
    (h : Φ.IsNonDilatingOn) : (phiScOn ℝ≥0 Φ hcharInj).IsNonDilatingOn := by
  refine MonoidOn.scOn_isNonDilatingOn Φ hcharInj
    (fun A => ((hpf A).choose_spec.divisorial).2)
    (fun A => isSharp_scT (hpf A).choose_spec (hpf A).choose_spec.divisorial)
    ?_ ?_ h
  · intro A b hb
    refine hprimWeak_of_isPrimaryElt (fun b' hb' => ?_) b hb
    exact isPrimaryElt_toSc (hpf A).choose_spec (hpf A).choose_spec.divisorial
      (suppElt_eq_singleton_toPrime (hpf A).choose_spec (hperf A)
        (hpf A).choose_spec.divisorial hb') hb'.1
  · intro A b hb x hx
    exact mprec_of_mprec_sc_primary_of_perfect (hpf A).choose_spec (hperf A)
      (hpf A).choose_spec.divisorial hb x hx

end ReflOn

/-! ## ★6. `Proposition 5.5, (iii)` の `𝒞^rlf` の側を組み上げる -/

section StdRlf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (iii) の `𝒞^rlf` の側(組み上げ)** ——
**`𝒞` が standard 型かつ not group-like で `Φ` が perf-factorial なら、
`𝒞^rlf` も standard 型**。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★`scModel_standardType` が要求していた 2 条が**両方とも消えた**:

| 条 | 出どころ |
|---|---|
| `hnd`(`Φ^rlf` が non-dilating) | ★`MonoidOn.scOn_isNonDilatingOn_of_perfFactorial`(本ファイル) |
| `hngl`(`𝒞^rlf` が not group-like) | ★在庫 `scModel_not_isOfGroupLikeType` ＋ 本ファイルの `isSharp_scT` / `toSc_injective` |

★★★残る仮引数は**原文の常備仮定**(`hiso` / `hfn` / `hcharInj` / `hint` / `hfsmD` /
`hdiv` / `htot` / `hconn` / `hfsmff` / `hpf` / `hnd` / `hngl`)と、
在庫補題が要求する `hperf` だけである。 -/
theorem scModel_standardType_of_perfFactorial (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := ℝ≥0) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT ℝ≥0 (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdiv : (scModel ℝ≥0 G hiso hfn hcharInj hint hfsmD).phi.IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D)
    (F' : FrobenioidCore (ModelData.modelPre
      (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)))
    (hfsmff : IsOfFSMFFType D)
    (hpf : ∀ A : D, IsPerfFactorial (Φ.val A))
    (hperf : ∀ A : D, IsPerfectMonoid (Φ.val A))
    (hnd : Φ.IsNonDilatingOn)
    (hngl : ¬ IsOfGroupLikeType P) :
    IsOfStandardType D (ScModelObj ℝ≥0 G hiso hfn hcharInj hint hfsmD)
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)) F' :=
  scModel_standardType G hiso hfn hcharInj hint hfsmD hdiv htot hconn F' hfsmff
    (MonoidOn.scOn_isNonDilatingOn_of_perfFactorial Φ hcharInj hpf hperf hnd)
    (scModel_not_isOfGroupLikeType G hiso hfn hcharInj hint hfsmD hdiv htot hconn
      (fun d => isSharp_scT (hpf d).choose_spec (hpf d).choose_spec.divisorial)
      (fun d => toSc_injective (hpf d).choose_spec (hpf d).choose_spec.divisorial)
      hngl)

end StdRlf

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
