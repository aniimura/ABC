/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38Kummer
import Mathlib.FieldTheory.KummerExtension
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★Kummer の `σ` は取れる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★訂正 —— `§9-994` の「mathlib に Kummer 対応は無い」は**誤り**だった

`§9-994` の `.needs` に「Kummer 対応そのものは mathlib に無い（2026-08-29 実測）」と
書いたが、★★**それは誤測定である**——見ていたのは `AlgebraicGeometry/` と
`NumberTheory/` だけで、**`Mathlib/FieldTheory/KummerExtension.lean` を見ていなかった**。

★mathlib には次がある（2026-08-29 に `#check` で確認）:

| 補題 | 内容 |
|---|---|
| `X_pow_sub_C_irreducible_of_prime` | `p` 素・`∀ b, b^p ≠ a` ⟹ `X^p − C a` は既約 |
| `autAdjoinRootXPowSubCEquiv` | `rootsOfUnity n K ≃* Gal(AdjoinRoot (X^n − C a) / K)` |
| ★`autAdjoinRootXPowSubCEquiv_root` | **`σ_η(root) = η • root`** |
| `X_pow_sub_C_splits_of_isPrimitiveRoot` | 分解体になること |

★★★**`autAdjoinRootXPowSubCEquiv_root` が求めていた `σ` そのもの**である。

## ★★★★★★★★★★★★★★★★★★★★これで取れること

    `l ∤ v_K(q)`  ⟹  `∀ b : K, b^l ≠ q`  ⟹  `X^l − C q` は既約
                  ⟹  各 `l` 乗根 `η` に対し `σ(π) = η·π` なる `σ` がある

★`§9-993` の `sigma_acts_as_alpha` と合わせると、`η` を原始的に取れば
**mod `l` 像が `α = (1 1 / 0 1)` を含む**。

## ★★これで `Theorem 3.8` に残るのは 3 点

| 段 | 状態 |
|---|---|
| `α` が像に入る | ★★**本ファイル ＋ `§9-993` ＋ `§9-994`** |
| `l`-巡回 ⟷ 安定直線の対応 | ☆残る（有限平坦群スキーム） |
| `Lemma 3.7` | ☆`Prop 3.4` 待ち |
| `torsionExt`（3・5 捩れ有理化 ⟹ 半安定） | ☆残る |
| 群論（`Lemma 3.1, (iv)`・`§9-992`） | ★済み |

★☆なお本ファイルは `AdjoinRoot (X^l − C q)` の上で述べている。
Tate 曲線の `L` へ移す段（`π` を `L` の中で取り、`Gal(L/K)` の元に持ち上げる）は
`Lemma32Tate.lean` の設定に合わせて別に要る。★**逸脱として明示する。**

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

open Polynomial

/-! ## ★★★★★★★★付値から「`l` 乗でない」へ -/

/-- ★★★★★★★★**`l ∤ v(a)` なら `a` は `K` で `l` 乗でない**。

★`b^l = a` なら `b ≠ 0` なので単元にでき、`l · v(b) = v(a)` から `l ∣ v(a)`。
★★`§9-994` の `not_lth_power_of_not_dvd_val`（分岐指数つきの版）の、
**基礎体の上での形**である。 -/
theorem not_lth_power_of_val {K : Type} [Field K] {l : ℕ} (hl : Nat.Prime l)
    (v : Kˣ →* Multiplicative ℤ) (a : Kˣ)
    (hnd : ¬ ((l : ℤ) ∣ Multiplicative.toAdd (v a))) :
    ∀ b : K, b ^ l ≠ (a : K) := by
  intro b hb
  have hb0 : b ≠ 0 := by
    intro h
    rw [h, zero_pow hl.ne_zero] at hb
    exact a.ne_zero hb.symm
  set u : Kˣ := Units.mk0 b hb0 with hu
  have hual : u ^ l = a := by
    ext
    push_cast
    simpa [hu] using hb
  have hv : v (u ^ l) = v a := by rw [hual]
  rw [map_pow] at hv
  refine hnd ⟨Multiplicative.toAdd (v u), ?_⟩
  have h := congrArg Multiplicative.toAdd hv
  simpa [mul_comm] using h.symm

/-! ## ★★★★★★★★★★★★★★★★★★★★Kummer の `σ` -/

/-- ★★★★★★★★★★★★★★★★★★★★**`σ(π) = η·π` なる `σ` がある**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `autAdjoinRootXPowSubCEquiv_root` そのものである
（`Mathlib/FieldTheory/KummerExtension.lean`）。
★★`§9-994` に「mathlib に Kummer 対応は無い」と書いたのは**誤測定**であった。 -/
theorem exists_sigma_smul_root {K : Type} [Field K] {l : ℕ} (hl : Nat.Prime l) [NeZero l]
    (hζ : (primitiveRoots l K).Nonempty) {a : K} (ha : ∀ b : K, b ^ l ≠ a)
    (η : rootsOfUnity l K) :
    ∃ σ : AdjoinRoot (X ^ l - C a) ≃ₐ[K] AdjoinRoot (X ^ l - C a),
      σ (AdjoinRoot.root (X ^ l - C a))
        = ((η : Kˣ) : K) • AdjoinRoot.root (X ^ l - C a) := by
  refine ⟨autAdjoinRootXPowSubCEquiv hζ (X_pow_sub_C_irreducible_of_prime hl ha) η, ?_⟩
  exact autAdjoinRootXPowSubCEquiv_root hζ _ η

/-- ★★★★★★★★★★★★★★★★★★★★★★**付値から直接**——
`l ∤ v_K(q)` なら各 `l` 乗根 `η` に対し `σ(π) = η·π` なる `σ` がある。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが原文 p.20 の「`the image of Galois in GL₂(𝔽_l)` contains the element `α`」の
**存在の側**である（行列表示は `§9-993`）。 -/
theorem exists_sigma_smul_root_of_val {K : Type} [Field K] {l : ℕ} (hl : Nat.Prime l)
    [NeZero l] (hζ : (primitiveRoots l K).Nonempty)
    (v : Kˣ →* Multiplicative ℤ) (q : Kˣ)
    (hnd : ¬ ((l : ℤ) ∣ Multiplicative.toAdd (v q)))
    (η : rootsOfUnity l K) :
    ∃ σ : AdjoinRoot (X ^ l - C (q : K)) ≃ₐ[K] AdjoinRoot (X ^ l - C (q : K)),
      σ (AdjoinRoot.root (X ^ l - C (q : K)))
        = ((η : Kˣ) : K) • AdjoinRoot.root (X ^ l - C (q : K)) :=
  exists_sigma_smul_root hl hζ (not_lth_power_of_val hl v q hnd) η

/-! ## ★★★★★★★★★★★★Kummer の根を単数として取る -/

open AdjoinRoot in
/-- ★★★★★★★★★★★★**Kummer の根は単数で、その `l` 乗は `q`**——★**無条件**（第 1178）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`q ≠ 0` なら `AdjoinRoot (Xˡ − C q)` の根 `π` は `0` でないので単数であり、
`πˡ = q`（`root_X_pow_sub_C_pow`）である。
★★★これが `tate_sigma_coord_alpha_of_base_root`（第 1177）が要求する
`π : Kˣ` と `πˡ = q` を**実際に取る**段である。 -/
theorem exists_kummer_root_unit {K : Type*} [Field K] {l : ℕ} (hl : 0 < l) {q : K}
    (hq : q ≠ 0) :
    ∃ π : (AdjoinRoot (X ^ l - C q))ˣ,
      (π : AdjoinRoot (X ^ l - C q)) = AdjoinRoot.root (X ^ l - C q) ∧
      (π : AdjoinRoot (X ^ l - C q)) ^ l = algebraMap K (AdjoinRoot (X ^ l - C q)) q := by
  have hpow : (AdjoinRoot.root (X ^ l - C q)) ^ l
      = algebraMap K (AdjoinRoot (X ^ l - C q)) q := root_X_pow_sub_C_pow l q
  have hqu : IsUnit (algebraMap K (AdjoinRoot (X ^ l - C q)) q) :=
    (isUnit_iff_ne_zero.2 hq).map (algebraMap K (AdjoinRoot (X ^ l - C q)))
  have hu : IsUnit (AdjoinRoot.root (X ^ l - C q)) := by
    obtain ⟨m, rfl⟩ : ∃ m, l = m + 1 := ⟨l - 1, by omega⟩
    have h1 : IsUnit ((AdjoinRoot.root (X ^ (m + 1) - C q)) ^ (m + 1)) := by
      rw [hpow]; exact hqu
    have h2 : (AdjoinRoot.root (X ^ (m + 1) - C q)) ^ (m + 1)
        = (AdjoinRoot.root (X ^ (m + 1) - C q)) ^ m
          * (AdjoinRoot.root (X ^ (m + 1) - C q)) := pow_succ _ m
    rw [h2] at h1
    exact isUnit_of_mul_isUnit_right h1
  refine ⟨hu.unit, hu.unit_spec, ?_⟩
  rw [hu.unit_spec]
  exact hpow

/-! ## ★★★★★★★★Kummer 拡大は体である -/

open AdjoinRoot in
/-- ★★★★★★★★★★**`q` が `l` 乗でなければ `AdjoinRoot (Xˡ − C q)` は体**——★**無条件**（第 1200）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`l` 素数で `q` が `l` 乗でなければ `Xˡ − C q` は既約
（mathlib の `X_pow_sub_C_irreducible_of_prime`）なので、`AdjoinRoot` は体になる。

★★★これで Kummer 拡大を **Tate 設定の体 `K`** として使えるようになる
——`TateSetup` は `[Field K]` を要求するからである。 -/
theorem irreducible_X_pow_sub_C_of_not_pow {K : Type*} [Field K] {l : ℕ} (hl : l.Prime)
    {q : K} (hq : ∀ b : K, b ^ l ≠ q) :
    Irreducible (X ^ l - C q) :=
  X_pow_sub_C_irreducible_of_prime hl hq

/-- ☆`Fact` の形（インスタンス解決に使う）。 -/
theorem fact_irreducible_X_pow_sub_C_of_not_pow {K : Type*} [Field K] {l : ℕ} (hl : l.Prime)
    {q : K} (hq : ∀ b : K, b ^ l ≠ q) :
    Fact (Irreducible (X ^ l - C q)) :=
  ⟨irreducible_X_pow_sub_C_of_not_pow hl hq⟩

/-! ## ★★★★★★★★★★★★★★★★Kummer の `σ` を単数の準同型として取る -/

open AdjoinRoot in
/-- ★★★★★★★★★★★★★★★★★★★★
**Kummer の `σ` を単数の準同型として取る**——★**無条件**（第 1211）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`l ∤ v(q)` なら `K = AdjoinRoot (Xˡ − C q)` の上に `σ(π) = ζ·π` なる
`K₀`-代数同型がある（第 994）。本定理はそれを

* `π : Kˣ`（`πˡ = q`、第 1178）
* `σ : Kˣ →* Kˣ`（`σ ζ = ζ`・`σ π = ζ·π`）

の形に**まとめる**——これが `tate_sigma_coord_alpha`（第 1174）が
受け取る形そのものである。

★★★これで `Skeleton/GenEll/AlphaBridge.lean` の
「Kummer の `σ` を Tate 設定に載せる」段の**単数の側**が済む。
☆残るのは `TateSetup` を `K` へ底変換する段である。 -/
theorem exists_units_sigma_kummer {K₀ : Type} [Field K₀] {l : ℕ} (hl : Nat.Prime l)
    (v : K₀ˣ →* Multiplicative ℤ) (q : K₀ˣ)
    (hnd : ¬ ((l : ℤ) ∣ Multiplicative.toAdd (v q)))
    (ζ₀ : K₀ˣ) (hζ₀ : IsPrimitiveRoot ((ζ₀ : K₀)) l) :
    ∃ (π : (AdjoinRoot (X ^ l - C (q : K₀)))ˣ)
      (σ : (AdjoinRoot (X ^ l - C (q : K₀)))ˣ →* (AdjoinRoot (X ^ l - C (q : K₀)))ˣ),
      π ^ l = Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) q ∧
      σ (Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) ζ₀)
        = Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) ζ₀ ∧
      σ π = Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) ζ₀ * π := by
  haveI : NeZero l := ⟨hl.ne_zero⟩
  have hζne : (primitiveRoots l K₀).Nonempty :=
    ⟨(ζ₀ : K₀), (mem_primitiveRoots hl.pos).2 hζ₀⟩
  have hζpow : ζ₀ ^ l = 1 := Units.ext (by simpa using hζ₀.pow_eq_one)
  let η : rootsOfUnity l K₀ := ⟨ζ₀, (mem_rootsOfUnity l ζ₀).2 hζpow⟩
  obtain ⟨τ, hτ⟩ := exists_sigma_smul_root_of_val hl hζne v q hnd η
  obtain ⟨π, hπ, hπl⟩ := exists_kummer_root_unit hl.pos q.ne_zero
  refine ⟨π, Units.map τ.toAlgHom.toRingHom.toMonoidHom, ?_, ?_, ?_⟩
  · exact Units.ext (by simpa using hπl)
  · exact Units.ext (by simpa using τ.commutes (ζ₀ : K₀))
  · refine Units.ext ?_
    show τ ((π : AdjoinRoot (X ^ l - C (q : K₀))))
      = algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) (ζ₀ : K₀)
        * (π : AdjoinRoot (X ^ l - C (q : K₀)))
    rw [hπ, hτ, Algebra.smul_def]

/-! ## ★★★★★★★★★★★★第 1274 —— 体自己同型として取り出す -/

/-- ★★★★★★★★★★★★★★★★
**Kummer の `σ` を体自己同型として取る**——★**無条件**（第 1274）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1211（`exists_units_sigma_kummer`）は単数群の準同型だけを返していたが、
`tatePhi_pointMap`（在庫）は**体の準同型 `σA : K →ₐ[R] K`** を要求する。
★そこで元の `AlgEquiv` も一緒に返す形にした——`σU ≔ Units.map τ` と取れば
`hσU : ↑(σU u) = σA ↑u` は `rfl` である。

★★★これで第 1273（`tate_point_unipotent`）に渡せる形が揃う。 -/
theorem exists_algEquiv_sigma_kummer {K₀ : Type} [Field K₀] {l : ℕ} (hl : Nat.Prime l)
    (v : K₀ˣ →* Multiplicative ℤ) (q : K₀ˣ)
    (hnd : ¬ ((l : ℤ) ∣ Multiplicative.toAdd (v q)))
    (ζ₀ : K₀ˣ) (hζ₀ : IsPrimitiveRoot ((ζ₀ : K₀)) l) :
    ∃ (π : (AdjoinRoot (X ^ l - C (q : K₀)))ˣ)
      (τ : AdjoinRoot (X ^ l - C (q : K₀)) ≃ₐ[K₀] AdjoinRoot (X ^ l - C (q : K₀))),
      π ^ l = Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) q ∧
      Units.map τ.toAlgHom.toRingHom.toMonoidHom π
        = Units.map (algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) : K₀ →* _) ζ₀ * π := by
  haveI : NeZero l := ⟨hl.ne_zero⟩
  have hζne : (primitiveRoots l K₀).Nonempty :=
    ⟨(ζ₀ : K₀), (mem_primitiveRoots hl.pos).2 hζ₀⟩
  have hζpow : ζ₀ ^ l = 1 := Units.ext (by simpa using hζ₀.pow_eq_one)
  let η : rootsOfUnity l K₀ := ⟨ζ₀, (mem_rootsOfUnity l ζ₀).2 hζpow⟩
  obtain ⟨τ, hτ⟩ := exists_sigma_smul_root_of_val hl hζne v q hnd η
  obtain ⟨π, hπ, hπl⟩ := exists_kummer_root_unit hl.pos q.ne_zero
  refine ⟨π, τ, Units.ext (by simpa using hπl), ?_⟩
  refine Units.ext ?_
  show τ ((π : AdjoinRoot (X ^ l - C (q : K₀))))
    = algebraMap K₀ (AdjoinRoot (X ^ l - C (q : K₀))) (ζ₀ : K₀)
      * (π : AdjoinRoot (X ^ l - C (q : K₀)))
  rw [hπ, hτ, Algebra.smul_def]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def irreducible_X_pow_sub_C_of_not_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(q が l 乗でなければ AdjoinRoot (Xˡ − C q) は体。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_algEquiv_sigma_kummer.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Kummer の σ を体自己同型として取る。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_units_sigma_kummer.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Kummer の σ を単数の準同型として取る。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_units_sigma_kummer.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_sigma_smul_root_of_val(第 994、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_sigma_smul_root_of_val") 1,
    .citation "[ABC3]" "exists_kummer_root_unit(第 1178、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_kummer_root_unit") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1211）**——`tate_sigma_coord_alpha`（第 1174）が" ++
       "受け取る形（`π : Kˣ`で `πˡ = q`、`σ : Kˣ →* Kˣ` で `σ ζ = ζ`・`σ π = ζ·π`）" ++
       "に**まとめた**ものである。" ++
       "☆残るのは `TateSetup` を `K` へ底変換する段である。") 3 ]

def exists_kummer_root_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Kummer の根は単数で、その l 乗は q。★無条件)",
    sectionId := "genell-thm-3-8" }

def not_lth_power_of_val.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l ∤ v(q) なら q は K で l 乗でない)",
    sectionId := "genell-thm-3-8" }

def exists_sigma_smul_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Kummer の σ——σ(π) = η·π)",
    sectionId := "genell-thm-3-8" }

def exists_sigma_smul_root_of_val.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(α が像に入る——存在の側)",
    sectionId := "genell-thm-3-8" }

def exists_sigma_smul_root_of_val.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "autAdjoinRootXPowSubCEquiv_root(σ_η(root) = η • root)"
      (.inMathlib "autAdjoinRootXPowSubCEquiv_root") 3,
    .citation "[mathlib]" "X_pow_sub_C_irreducible_of_prime(X^p − C a の既約性)"
      (.inMathlib "X_pow_sub_C_irreducible_of_prime") 2,
    .implicitStep
      ("★★★★★★★★★★訂正(2026-08-29): §9-994 の .needs に" ++
       "「Kummer 対応そのものは mathlib に無い」と書いたのは**誤測定**である" ++
       "——見ていたのは AlgebraicGeometry/ と NumberTheory/ だけで、" ++
       "**Mathlib/FieldTheory/KummerExtension.lean を見ていなかった**。" ++
       "★autAdjoinRootXPowSubCEquiv_root が求めていた σ そのものであった") 7,
    .implicitStep
      ("★★★逸脱: 本ファイルは AdjoinRoot (X^l − C q) の上で述べている。" ++
       "Tate 曲線の L へ移す段(π を L の中で取り、Gal(L/K) の元に持ち上げる)は" ++
       "Lemma32Tate.lean の設定に合わせて別に要る") 5,
    .implicitStep
      ("★★Theorem 3.8 の残高(2026-08-29): " ++
       "(1) l-巡回 ⟷ 安定直線の対応(有限平坦群スキーム)、" ++
       "(2) Lemma 3.7(Prop 3.4 待ち)、(3) torsionExt。" ++
       "★『α が像に入る』は §9-993(行列表示)＋§9-994(付値の障害)＋本ファイル(存在)で取れた。" ++
       "★★群論(Lemma 3.1, (iv)・§9-992)と Galois 表現の構成(galRep)も済んでいる") 7 ]

end ABC3.Found.GenEll
