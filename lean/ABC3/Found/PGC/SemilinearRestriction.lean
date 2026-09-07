import ABC3.Found.PGC.LubinTateTowerTransport

/-!
# 半線型な輸送の「制限」——(ii) の 1 と 3、および Lubin-Tate 作用の同変性

直前の節点(`LubinTateTowerTransport.lean`)は [pGC] Proposition 1.1 の残りを
(i)(ii)(iii) の 3 段に分け、**(i) を完成**させたうえで、(ii) に足りないものを
3 つ名指しした:

> **1.** `Φ` を `adjoinIntegers K x → adjoinIntegers K (Φ x)` に制限する段
>   ((a) ノルム保存、(b) 制限の連続性)。
>   ★これが木が「cross-point instance bridging」と呼んで避けてきた段である。
> **2.** `reciprocityUnits` の捩れ点の上での spec。
> **3.** 生成元列の付け替え(任意の生成元に移す段)。

本ファイルは **1 を完全に埋め**、その上で

* ★★**`Φ ([a]_f(x)) = [σ a]_{f^σ}(Φ x)`**(`semilinear_lubinTateActionAtTorsionPoint`)
* ★★**`ρ_{f^σ, n}(Φ τ Φ^{-1}) = σ(ρ_{f,n}(τ))`**(`reciprocityMap_semilinear_conj`)

を証明した。★**後者が申し送りの (ii) の有限段(レベル `n`)そのもの**である。
★2 と 3 は**この有限段の中では要らなかった**——理由は下の「段取りとの差分」に書く。

★★**`ArtinUnitEquivariance` は出ていない。**出ていないものを出たと書かない。
何が残っているかは末尾の「残っているものを名指しで」に書いた。

## ★★何がこのファイルで出たか(★測定結果)

| 段 | 宣言 | 状態 |
|---|---|---|
| 核 A | `spectralNorm_semilinear_ringEquiv` | ★★**証明した** —— 半線型な環同型はスペクトルノルムを保つ |
| 核 B | `restrictIntermediateFieldEquiv` / `restrictSubringEquiv` | ★**証明した**(純代数) |
| 核 C | `conjSemilinearAlgEquiv` | ★**証明した** —— `Φ τ Φ^{-1}` は `k`-代数同型 |
| 核 D | `map_iteratedLubinTatePsi` | ★**証明した** —— `ψ_n` の係数ひねりの自然性 |
| 1-(a) | `norm_semilinear_closureEquiv` | ★★**証明した**(核 A に代入するだけ) |
| 1-(b) | `continuous_adjoinIntegersSemilinear` | ★★**証明した**(等長 ⇒ 連続) |
| 1 | `adjoinIntegersSemilinear` | ★★★**証明した** —— これが「cross-point instance bridging」 |
| (ii-点) | `semilinear_lubinTateActionAtTorsionPoint` | ★★★★**証明した** |
| (ii-n) | `reciprocityMap_semilinear_conj` | ★★★★★**証明した** —— 申し送りの (ii) の有限段 |
| 橋 | `absGalConjCME_eq_conjSemilinearAlgEquiv` | ★**証明した(`rfl`)** —— `Ψ_g` = 抽象核 C |
| 副産物 | `algEquiv_lubinTateActionAtTorsionPoint_cross` | ★★★★**証明した** —— **点をまたぐ** Galois 同変性 |
| (iii) | 素元非依存 | ★**要らなかった**(有限段では。下の「段取りとの差分」参照) |
| 壁 | `ArtinUnitEquivariance` | ★★**出ていない** |

## ★★★「cross-point instance bridging」はどう越えたか(★本持ち場の重心)

申し送りは「(a) `Φ` のノルム保存(★木の `norm_algEquiv_eq` は `≃ₐ[K.carrier]` 専用で
半線型版は無い)、(b) 制限の連続性」を壁として名指ししていた。
★★**測ったところ、両方が 1 本の抽象核から出た。**

* (a) は **mathlib の `spectralNorm_unique_field_norm_ext`**(完備な基点の上では
  基点のノルムを延長する絶対値は spectralNorm しかない)に、
  「`spectralNorm k Ω (Φ ·)` は絶対値である」を代入するだけで出る。
  ★`Φ` の `k`-線型性を**一切使わない**。使うのは
  「`Φ` は基点の上でノルムを保つ半線型」だけである。
  ★基点の上のノルム保存(`∀ c : K.carrier, ‖σ c‖ = ‖c‖`)は
  `σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier` なので mathlib の
  `spectralNorm_eq_of_equiv` が 1 行で与える(K.carrier のノルムは
  `spectralNorm ℚ_[p] K.carrier` そのものだから)。
* (b) は (a) の**系**である:ノルムを保つ加法準同型は等長
  (`AddMonoidHomClass.isometry_of_norm`)、等長なら連続。
  ★`LubinTateActionEquivariance.lean` の
  `continuous_algEquivRestrictSelf` は「有限次元だから線型写像は連続」
  (`LinearMap.continuous_of_finiteDimensional`)を使うが、
  ★**半線型ではその道は塞がっている**(`Φ` は `K.carrier`-線型でない)。
  ★等長性の道は**基礎体を一切見ない**ので半線型でもそのまま通る。

★★**したがって「cross-point instance bridging」は避けられなかったが、越えられた。**
★★**そして越えるのに要ったのは「等長性」だけで、`Φ x ∈ K⟮x⟯` のような
留まる性質は 1 度も使っていない**(実際 `Φ x` は `f^σ` の捩れ点であって
`f` の捩れ点ではないので、留まらない)。

## ★★段取りとの差分(★「思ったより安い」も同じ価値)

* ★★**2(`reciprocityUnits` の spec)と 3(生成元列の付け替え)は、
  有限段 `ρ_{f,n}` の同変性には要らなかった。**
  理由:木の `reciprocityMap`(レベル `n`、点 `x` を固定した版)は
  **定義そのものが「`x` の上での記述」**(`reciprocityMap_spec`)であり、
  一意性(`existsUnique_unitActionQuotient_eq_algEquiv`)も既に在庫にある。
  ★2 と 3 が要るのは**極限**(`reciprocityUnits`、`psiGenSeq` を使う版)へ
  上げるときであって、有限段では要らない。
* ★★**`LubinTateEndoTwisted` は要らなかった**(直前の波と同じ理由:
  `f` と `f^σ` は同じ環の上なので untwisted の一意性で足りる)。
  ★本持ち場では `powerSeries_uniqueness` すら直接は使っていない
  ——`map_LubinTateAction`(直前の波)に代入するだけで済んだ。
* ★★**素元非依存性(iii)は 1 度も使わなかった。**★理由を測った:
  有限段の (ii) は `ρ_{f,n}`(`f` の塔の上)と `ρ_{f^σ,n}`(`f^σ` の塔の上)を
  **それぞれ自分の塔の上で**比べる形なので、2 つの塔が一致する必要が無い。
  ★**(iii) が要るのは「同じ `ρ` である」と言いたいとき**——つまり
  `ArtinUnitEquivariance` の `E` を作る段である。
  ★直前の波の測定(`ρ_f ≠ ρ_{f^{σ_g}}`、一致するのは慣性部分＝`p^n` 捩れの上だけ)は
  そのまま生きている。
* ★★**副産物: 「点をまたぐ Galois 同変性」(§10)が落ちた。**
  `σ := AlgEquiv.refl` を入れるだけで `τ(a·x) = a·(τ x)` が
  **`adjoinIntegers K (τ x)` の中の本物の作用**として出る。
  ★木は 20 波以上これを避けて `lubinTateActionAtAlgEquivPoint`
  (`x` の座標系での代用品)で回していた。★**`τ x ∈ K⟮x⟯` を使っていないので
  留まらない場合にもそのまま成り立つ。**
* ★**`ψ_n` の自然性(`map_iteratedLubinTatePsi`)が新しく要った。**
  `ψ_n` は `Exists.choose` を経由して定義されている(`iteratedLubinTateStepR`)ので
  自然性が定義から出ない。★**`D_n = D_{n-1}·ψ_n`(在庫)+ `D_n` の自然性(在庫)+
  多項式環が整域であること、の 3 つで割り算して出した**(20 行)。
* ★**#59 には 1 度も当たっていない。** 中間体は `K⟮x⟯` の 1 層だけで、
  `restrictScalars` も `restrictNormalHom` も書いていない(#59 の定型 (d))。

## ★★退化の自己検査

* ★**`Φ x` は `K⟮x⟯` に留まらない。** 本ファイルはどこでも留まると仮定していない
  ——`adjoinIntegersSemilinear` の終域は `adjoinIntegers K (Φ x)` という**別の環**である。
* ★**`Φ` に `K.carrier`-線型性を仮定していない。** 仮定すると `σ = id` になり主張が空になる。
* ★**基点は `ℚ_[p]`。** `σ` は `≃ₐ[ℚ_[p]]` であって `≃ₐ[K.carrier]` ではない。
* ★**`reciprocityMap_semilinear_conj` は `u` の代表元を明示的に取る形**にした
  (`Units.map` を書かない)。★これは #205(`⇑↑Φ` と `⇑Φ` の食い違い)を
  statement の段階で避けるためである。

## ★★★残っているものを名指しで(★次の節点はここから持てる)

1. ★**極限への持ち上げ**:`reciprocityMap_semilinear_conj`(レベル `n`)を
   `reciprocityUnits`(`LubinTateClosureTopology.lean`)へ上げる段。
   ★ここで初めて申し送りの **2**(`reciprocityUnits` の捩れ点上の spec)と
   **3**(`psiGenSeq` の付け替え)が要る。
   ★本ファイルの `reciprocityMap_semilinear_conj` は**任意の**
   `x`(`ψ_n` の根)について成り立つので、3 の「任意生成元に移す段」は
   もう **statement の側で済んでいる**。
2. ★**(iii) 素元非依存**:`Φ x` は `f^σ` の捩れ点なので、塔は `π` から `σ(π)` へ
   動く。★直前の波の測定どおり `ρ_f ≠ ρ_{f^σ}` であって、
   一致するのは慣性部分(`p^n` 捩れが乗る部分)だけ。
   ★**本ファイルはそこに触れていない**(有限段の `ρ_{f,n}` は
   `f` と `f^σ` の**別々の**塔の上で比べているので、素元非依存は要らなかった)。
3. ★**`Gal(L^{ab}/L) ≅ 𝒪_L^× × Ẑ` の分解 `E` を作る段**。
   `ArtinUnitEquivariance` は `E` の**存在**を要求するので、
   `reciprocityUnits` から `abelianClosure` の Galois 群への同型
   (`nonempty_abelianGalContinuousEquivUnitsZHat` が与えるのは
   ★**同型の存在だけで、それが `ρ` から来ていることは言っていない**)を
   `ρ` 由来のものに取り替える段が要る。★**ここが次の壁である。**
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

open scoped NormedField Valued Classical

/-! ## 1. 抽象核 A —— 半線型な環同型はスペクトルノルムを保つ

★この節に分岐・Lubin-Tate・Galois の語彙は 1 つも出てこない。
一般の完備非アルキメデス体 `k` とその代数拡大 `Ω` の話である。 -/

section AbstractCoreA

variable {k Ω : Type*} [NontriviallyNormedField k] [IsUltrametricDist k] [CompleteSpace k]
  [Field Ω] [Algebra k Ω] [Algebra.IsAlgebraic k Ω]

/-- ★**抽象核 A1** —— 環同型 `Φ` に沿ったスペクトルノルムの引き戻しは絶対値。

★`spectralMulAlgNorm`(mathlib)が乗法性・劣加法性を与え、
`eq_zero_of_map_spectralNorm_eq_zero` が非退化性を与える。
`Φ` が全単射であることは非退化性にだけ使う。 -/
noncomputable def pullbackSpectralAbs (Φ : Ω ≃+* Ω) : AbsoluteValue Ω ℝ where
  toFun y := spectralNorm k Ω (Φ y)
  map_mul' y z := by
    rw [map_mul, ← spectralMulAlgNorm_def (K := k) (L := Ω),
      ← spectralMulAlgNorm_def (K := k) (L := Ω) (Φ y),
      ← spectralMulAlgNorm_def (K := k) (L := Ω) (Φ z), map_mul]
  nonneg' y := spectralNorm_nonneg _
  eq_zero' y := by
    refine ⟨fun h => ?_, fun h => ?_⟩
    · have h0 := eq_zero_of_map_spectralNorm_eq_zero (K := k) (L := Ω) h
        (Algebra.IsAlgebraic.isAlgebraic _)
      exact Φ.injective (a₁ := y) (a₂ := 0) (by rw [h0, map_zero])
    · subst h
      show spectralNorm k Ω (Φ 0) = 0
      rw [map_zero, spectralNorm_zero]
  add_le' y z := by
    show spectralNorm k Ω (Φ (y + z)) ≤ spectralNorm k Ω (Φ y) + spectralNorm k Ω (Φ z)
    rw [map_add, ← spectralMulAlgNorm_def (K := k) (L := Ω),
      ← spectralMulAlgNorm_def (K := k) (L := Ω) (Φ y),
      ← spectralMulAlgNorm_def (K := k) (L := Ω) (Φ z)]
    exact map_add_le_add _ _ _

/-- ★★★★**抽象核 A —— 基点の上でノルムを保つ半線型な環同型は
スペクトルノルムを保つ**:`‖Φ y‖ = ‖y‖`。

★★**`Φ` に `k`-線型性を仮定していない**(仮定すると本ファイルの主張が空になる)。
仮定しているのは
* `Φ` は `σ` について半線型(`Φ ∘ ι = ι ∘ σ`)、
* `σ` は `k` のノルムを保つ、
の 2 つだけである。

★これが「cross-point instance bridging」の (a) を潰す核である。
木の `norm_algEquiv_eq`(`AdjoinIntegers.lean`)は `≃ₐ[K.carrier]` 専用で、
最小多項式の根の集合を経由するので半線型には効かない。
★こちらは **mathlib の `spectralNorm_unique_field_norm_ext`(一意性)** に
乗せるだけなので、線型性をどこにも使わない。 -/
theorem spectralNorm_semilinear_ringEquiv (Φ : Ω ≃+* Ω) (σ : k → k)
    (hσ : ∀ c : k, ‖σ c‖ = ‖c‖)
    (hΦ : ∀ c : k, Φ (algebraMap k Ω c) = algebraMap k Ω (σ c)) (y : Ω) :
    spectralNorm k Ω (Φ y) = spectralNorm k Ω y :=
  spectralNorm_unique_field_norm_ext (f := pullbackSpectralAbs Φ)
    (fun c => by
      show spectralNorm k Ω (Φ (algebraMap k Ω c)) = ‖c‖
      rw [hΦ, spectralNorm_extends, hσ]) y

end AbstractCoreA

/-! ## 2. 抽象核 B —— 制限(中間体・部分環)

★この節は純代数である(位相もノルムも出てこない)。 -/

section AbstractCoreB

/-- ★**抽象核 B1** —— 像が一致する 2 つの中間体の間への環同型の制限。

★`E` と `F` は**別の**中間体でよい(`Φ` が `E` を自分自身に写すことを仮定しない)。
これが半線型の場合に必要な形である。 -/
def restrictIntermediateFieldEquiv {Ω : Type*} [Field Ω] {k : Type*} [Field k] [Algebra k Ω]
    (Φ : Ω ≃+* Ω) (E F : IntermediateField k Ω)
    (h : Φ '' (E : Set Ω) = (F : Set Ω)) : (E : Type _) ≃+* (F : Type _) where
  toFun y := ⟨Φ (y : Ω), by
    show Φ (y : Ω) ∈ (F : Set Ω)
    rw [← h]
    exact ⟨(y : Ω), y.2, rfl⟩⟩
  invFun z := ⟨Φ.symm (z : Ω), by
    have hz : (z : Ω) ∈ Φ '' (E : Set Ω) := by rw [h]; exact z.2
    obtain ⟨y, hy, hyz⟩ := hz
    show Φ.symm (z : Ω) ∈ (E : Set Ω)
    rw [← hyz, Φ.symm_apply_apply]
    exact hy⟩
  left_inv y := Subtype.ext (Φ.symm_apply_apply (y : Ω))
  right_inv z := Subtype.ext (Φ.apply_symm_apply (z : Ω))
  map_mul' a b := Subtype.ext (by push_cast; rw [map_mul])
  map_add' a b := Subtype.ext (by push_cast; rw [map_add])

@[simp] theorem coe_restrictIntermediateFieldEquiv {Ω : Type*} [Field Ω] {k : Type*} [Field k]
    [Algebra k Ω] (Φ : Ω ≃+* Ω) (E F : IntermediateField k Ω)
    (h : Φ '' (E : Set Ω) = (F : Set Ω)) (y : E) :
    ((restrictIntermediateFieldEquiv Φ E F h y : F) : Ω) = Φ (y : Ω) := rfl

/-- ★**抽象核 B2** —— 部分環を部分環へ写す環同型の制限。 -/
def restrictSubringEquiv {A B : Type*} [Ring A] [Ring B] (e : A ≃+* B)
    (SA : Subring A) (SB : Subring B) (h : ∀ a : A, a ∈ SA ↔ e a ∈ SB) :
    (SA : Type _) ≃+* (SB : Type _) where
  toFun z := ⟨e (z : A), (h (z : A)).mp z.2⟩
  invFun w := ⟨e.symm (w : B), (h _).mpr (by rw [e.apply_symm_apply]; exact w.2)⟩
  left_inv z := Subtype.ext (e.symm_apply_apply (z : A))
  right_inv w := Subtype.ext (e.apply_symm_apply (w : B))
  map_mul' a b := Subtype.ext (by push_cast; rw [map_mul])
  map_add' a b := Subtype.ext (by push_cast; rw [map_add])

@[simp] theorem coe_restrictSubringEquiv {A B : Type*} [Ring A] [Ring B] (e : A ≃+* B)
    (SA : Subring A) (SB : Subring B) (h : ∀ a : A, a ∈ SA ↔ e a ∈ SB) (z : SA) :
    ((restrictSubringEquiv e SA SB h z : SB) : B) = e (z : A) := rfl

end AbstractCoreB

/-! ## 3. 抽象核 C —— 半線型な共役

★★原典 Corollary 4.9 の証明の 1 行 `σ([θ](α)) = [θ]^{(j)}[xπ_j](α)` は、
「`σ` を `[θ]` で共役する」という操作である。その**型の段**をここで切り出す。 -/

section AbstractCoreC

variable {k Ω : Type*} [CommRing k] [CommRing Ω] [Algebra k Ω]

/-- ★★**抽象核 C** —— 半線型な環同型による共役は `k`-代数同型を `k`-代数同型へ移す。

`Φ τ Φ^{-1}` は、`Φ` が `k` を(要素ごとにではなく)**集合として**動かすだけなので、
`k` を要素ごとに固定する。★これが「`Φ τ Φ^{-1}` が `Γ_K` に属する」ことの中身である。 -/
def conjSemilinearAlgEquiv (Φ : Ω ≃+* Ω) (σ : k ≃+* k)
    (hΦ : ∀ c : k, Φ (algebraMap k Ω c) = algebraMap k Ω (σ c)) (τ : Ω ≃ₐ[k] Ω) : Ω ≃ₐ[k] Ω :=
  AlgEquiv.ofRingEquiv (f := Φ.symm.trans (τ.toRingEquiv.trans Φ)) (fun c => by
    have hsymm : Φ.symm (algebraMap k Ω c) = algebraMap k Ω (σ.symm c) := by
      apply Φ.injective
      rw [Φ.apply_symm_apply, hΦ, σ.apply_symm_apply]
    show Φ (τ (Φ.symm (algebraMap k Ω c))) = algebraMap k Ω c
    rw [hsymm, τ.commutes, hΦ, σ.apply_symm_apply])

@[simp] theorem conjSemilinearAlgEquiv_apply (Φ : Ω ≃+* Ω) (σ : k ≃+* k)
    (hΦ : ∀ c : k, Φ (algebraMap k Ω c) = algebraMap k Ω (σ c)) (τ : Ω ≃ₐ[k] Ω) (y : Ω) :
    conjSemilinearAlgEquiv Φ σ hΦ τ y = Φ (τ (Φ.symm y)) := rfl

end AbstractCoreC

/-! ## 4. 抽象核 D —— `ψ_n` の係数ひねりの自然性

★★`ψ_n` は `iteratedLubinTateStepR`(`Exists.choose` で選ばれる)の
Weierstrass distinguished 部分として定義されているので、
**自然性が定義から出ない**。★`D_n = D_{n-1}·ψ_n`(在庫)と
`D_n` の自然性(在庫)から**割り算で**出す。

★この節は一般の完備局所整域 `A` の話で、体・付値・Galois の語彙は出てこない。 -/

section AbstractCoreD

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp]
    [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ}

/-- ★★**抽象核 D —— `ψ_n^{f^φ} = (ψ_n^f)^φ`**。

★在庫 `map_iteratedLubinTateDistinguished`(`D_n` の自然性、`ArtinEquivarianceProof.lean`)と
`iteratedLubinTateDistinguished_eq_mul_psi`(`D_n = D_{n-1}·ψ_n`、`LubinTateActionPsi.lean`)を
`Polynomial A` が整域であることで割るだけ。 -/
theorem map_iteratedLubinTatePsi
    (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff) (φ : A ≃+* A)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    iteratedLubinTatePsi hq (maximalIdeal_eq_span_map φ hπmax)
        (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
        (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
        (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n hn
      = Polynomial.map (φ : A →+* A)
          (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) := by
  have hD : ∀ k : ℕ,
      iteratedLubinTateDistinguished hq (maximalIdeal_eq_span_map φ hπmax)
          (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
          (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
          (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) k
        = Polynomial.map (φ : A →+* A)
            (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf k) := fun k =>
    map_iteratedLubinTateDistinguished hq φ hπmax hπne0 f hf0 hf1 hf k
  have hmul' := iteratedLubinTateDistinguished_eq_mul_psi hq (maximalIdeal_eq_span_map φ hπmax)
    (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
    (PowerSeries.map (φ : A →+* A) f) (coeff_zero_map_twist _ hf0)
    (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n hn
  have hmul := iteratedLubinTateDistinguished_eq_mul_psi hq hπmax hπne0 f hf0 hf1 hf n hn
  rw [hD n, hD (n - 1), hmul, Polynomial.map_mul] at hmul'
  have hne : Polynomial.map (φ : A →+* A)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)) ≠ 0 :=
    Polynomial.Monic.ne_zero
      ((isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1
        hf (n - 1)).monic.map _)
  exact mul_left_cancel₀ hne hmul'.symm

end AbstractCoreD

variable {p : ℕ} [Fact p.Prime]

/-! ## 5. ★★★★具体層 1 —— 「cross-point instance bridging」

★★申し送りの **1** がここで埋まる。 -/

/-- ★★**(1-a) 半線型な `Φ` は `K̄` のノルムを保つ**。

★抽象核 A に代入するだけ。基点の上のノルム保存は
`σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier` なので mathlib の `spectralNorm_eq_of_equiv`
(`K.carrier` のノルムは `spectralNorm ℚ_[p] K.carrier` そのもの)が与える。

★★木の `norm_algEquiv_eq` は `≃ₐ[K.carrier]` 専用で、これは**その半線型版**である。 -/
theorem norm_semilinear_closureEquiv (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (y : K.closure) :
    ‖Φ y‖ = ‖y‖ :=
  spectralNorm_semilinear_ringEquiv Φ σ (fun c => (spectralNorm_eq_of_equiv σ c).symm) hΦ y

/-- ★`Φ '' K⟮x⟯ = K⟮Φ x⟯` —— 直前の波の `image_adjoin_of_semilinear_equiv` に
`S := {x}` を入れるだけ。 -/
theorem image_adjoin_singleton_semilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) :
    Φ '' (IntermediateField.adjoin K.carrier ({x} : Set K.closure) : Set K.closure)
      = (IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure) : Set K.closure) := by
  rw [image_adjoin_of_semilinear_equiv Φ σ.toRingEquiv hΦ ({x} : Set K.closure),
    Set.image_singleton]

/-- ★★**`Φ` の `K⟮x⟯ → K⟮Φ x⟯` への制限**。

★★**終域は別の中間体である**(`Φ x` は `f^σ` の捩れ点なので `K⟮x⟯` に留まらない)。 -/
noncomputable def adjoinFieldSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) :
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure) : Type _) ≃+*
      (IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure) : Type _) :=
  restrictIntermediateFieldEquiv Φ _ _ (image_adjoin_singleton_semilinear K σ Φ hΦ x)

theorem coe_adjoinFieldSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure)
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ((adjoinFieldSemilinear K σ Φ hΦ x y :
        IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure)
      = Φ (y : K.closure) := rfl

/-- ★制限もノルムを保つ(部分体のノルムは環境のノルムだから)。 -/
theorem norm_adjoinFieldSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure)
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖adjoinFieldSemilinear K σ Φ hΦ x y‖ = ‖y‖ :=
  norm_semilinear_closureEquiv K σ Φ hΦ (y : K.closure)

/-- ★★**(1-b) の核** —— 制限は等長。 -/
theorem isometry_adjoinFieldSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) :
    Isometry (adjoinFieldSemilinear K σ Φ hΦ x) :=
  AddMonoidHomClass.isometry_of_norm _ (norm_adjoinFieldSemilinear K σ Φ hΦ x)

/-- ★★★★**申し送りの 1 そのもの** —— `Φ` の
`adjoinIntegers K x → adjoinIntegers K (Φ x)` への制限。

★★これが木が「cross-point instance bridging」と呼んで避けてきた段である。
半線型では避けられない(避けると `σ = id` になって主張が空になる)。
★越えるのに要ったのは**ノルム保存だけ**——`adjoinIntegers` は
「ノルム `≤ 1`」で定義されているので、ノルムを保てば自動的に対応する。 -/
noncomputable def adjoinIntegersSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) :
    (adjoinIntegers K x : Type _) ≃+* (adjoinIntegers K (Φ x) : Type _) :=
  restrictSubringEquiv (adjoinFieldSemilinear K σ Φ hΦ x) (adjoinIntegers K x)
    (adjoinIntegers K (Φ x))
    (fun a => by
      show ‖a‖ ≤ 1 ↔ ‖adjoinFieldSemilinear K σ Φ hΦ x a‖ ≤ 1
      rw [norm_adjoinFieldSemilinear])

theorem coe_adjoinIntegersSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) (z : adjoinIntegers K x) :
    (((adjoinIntegersSemilinear K σ Φ hΦ x z : adjoinIntegers K (Φ x)) :
        IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure)
      = Φ ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) := rfl

/-- ★★**(1-b) 制限は連続** —— 等長性の**系**である。

★`LubinTateActionEquivariance.lean` の `continuous_algEquivRestrictSelf` は
`LinearMap.continuous_of_finiteDimensional`(有限次元の線型写像は連続)を使うが、
★**半線型ではその道は塞がっている**。等長性の道は基礎体を見ないので通る。 -/
theorem continuous_adjoinIntegersSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) :
    Continuous (adjoinIntegersSemilinear K σ Φ hΦ x) := by
  apply Continuous.subtype_mk
  exact (isometry_adjoinFieldSemilinear K σ Φ hΦ x).continuous.comp continuous_subtype_val

/-- ★制限は `integerRingEquiv σ` について半線型。 -/
theorem semilinear_adjoinIntegersSemilinear (K : PAdicLocalField p)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a)) (x : K.closure) (c : 𝒪[K.carrier]) :
    adjoinIntegersSemilinear K σ Φ hΦ x (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) c)
      = algebraMap 𝒪[K.carrier] (adjoinIntegers K (Φ x)) (integerRingEquiv σ c) := by
  apply Subtype.ext
  apply Subtype.ext
  show Φ (algebraMap K.carrier K.closure (c : K.carrier))
      = algebraMap K.carrier K.closure (σ (c : K.carrier))
  exact hΦ (c : K.carrier)

/-! ## 6. 具体層 2 —— `ψ_n` の根の輸送

★直前の波の `map_mem_iteratedLubinTateTorsionPoints`(`Λ_n` の版)の `ψ_n` 版。
★**原始的な捩れ点であることが `Φ` で保たれる**——これが
`reciprocityMap` を `Φ x` の上で使うために要る。 -/

open scoped Classical in
/-- ★★**`Φ` は `ψ_n^f` の根を `ψ_n^{f^φ}` の根に写す**。

★証明の骨は `map_mem_iteratedLubinTateTorsionPoints` と同じで、
`D_n` の自然性の代わりに抽象核 D(`ψ_n` の自然性)を使うだけ。 -/
theorem map_mem_iteratedLubinTatePsiTorsionPoints (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (φ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier])
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (Φ : K.closure →+* K.closure)
    (hcompat : ∀ a : 𝒪[K.carrier], Φ (algebraMap 𝒪[K.carrier] K.closure a)
      = algebraMap 𝒪[K.carrier] K.closure (φ a))
    {x : K.closure}
    (hx : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    Φ x ∈ iteratedLubinTatePsiTorsionPoints K hq (maximalIdeal_eq_span_map φ hπmax)
      (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
      (PowerSeries.map (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) (coeff_zero_map_twist _ hf0)
      (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n hn := by
  rw [iteratedLubinTatePsiTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hx ⊢
  obtain ⟨-, hroot⟩ := hx
  refine ⟨?_, ?_⟩
  · refine Polynomial.Monic.ne_zero ?_
    exact ((isDistinguishedAt_iteratedLubinTatePsi hq (maximalIdeal_eq_span_map φ hπmax)
      (fun h => hπne0 (φ.injective (by rw [h, map_zero])))
      (PowerSeries.map (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) f) (coeff_zero_map_twist _ hf0)
      (coeff_one_map_twist _ hf1) (map_residue_map_twist φ hf) n hn).monic).map _
  · rw [Polynomial.IsRoot.def, Polynomial.eval_map] at hroot ⊢
    rw [map_iteratedLubinTatePsi hq φ hπmax hπne0 f hf0 hf1 hf n hn]
    exact eval₂_map_twist_eq_zero (φ : 𝒪[K.carrier] →+* 𝒪[K.carrier]) Φ
      (algebraMap 𝒪[K.carrier] K.closure) hcompat _ hroot

/-! ## 7. ★★★★払い出し A —— Lubin-Tate 作用の半線型同変性 -/

def semilinear_lubinTateActionAtTorsionPoint.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 8, item := "Lemma 4.6", sectionId := "lemma-4-6" }

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`Φ([a]_f(x)) = [σ a]_{f^σ}(Φ x)`** —— 半線型な `Φ` は
Lubin-Tate の `𝒪_K` 作用を「係数をひねった作用」に運ぶ。

原文 (Yoshida08 p.8):
> Lemma 4.6. Let f, f′ ∈ O[scr]_L[X] be as above with linear coefficients π, π′, respectively. If θ ∈ Θ^L,×_π,π′ (see Corollary 3.7(ii)), then for all m ≥ 1, it gives an isomorphism [θ] = [θ]_f,f′ : µ_f,m → µ_f′,m of O[scr]-modules, and L^m_f = L^m_f′.

★★**原典と同じ主張ではない**(逸脱として記録する)。原典は `[θ]`(`𝒪`-**線型**)による
輸送で `µ_{f,m} ≅ µ_{f′,m}` を出すが、ここでは `σ`(体の**半線型**自己同型)による
輸送で「`𝒪`-作用が係数のひねりと可換」を出す。★役割(「捩れ点の加群構造が
`f` の取り替えで運ばれる」段)が同じなので `.src` を Lemma 4.6 に付けた。

★証明は直前の波の抽象核 `aeval_powerSeries_comm_twist` に
本ファイルの `adjoinIntegersSemilinear` を代入し、冪級数の層
`map_LubinTateAction`(`([a]_f)^φ = [φa]_{f^φ}`)を使うだけである。 -/
theorem semilinear_lubinTateActionAtTorsionPoint (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a))
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (a : 𝒪[K.carrier]) :
    Φ ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)
      = ((lubinTateActionAtTorsionPoint K hq
            (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
            (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
            (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
              𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
            (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
            (map_residue_map_twist (integerRingEquiv σ) hf) n (Φ x)
            (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
              f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
              (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hx)
            (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x))
            (integerRingEquiv σ a) :
          IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  haveI := completeSpace_adjoinIntegers K (Φ x)
  haveI := isLinearTopology_adjoinIntegers K (Φ x)
  haveI := continuousSMul_adjoinIntegers K (Φ x)
  have hkey := aeval_powerSeries_comm_twist
    ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) : 𝒪[K.carrier] →+* 𝒪[K.carrier])
    ((adjoinIntegersSemilinear K σ Φ hΦ x : (adjoinIntegers K x : Type _) ≃+*
        (adjoinIntegers K (Φ x) : Type _)) :
      (adjoinIntegers K x : Type _) →+* (adjoinIntegers K (Φ x) : Type _))
    (continuous_adjoinIntegersSemilinear K σ Φ hΦ x)
    (semilinear_adjoinIntegersSemilinear K σ Φ hΦ x)
    (LubinTateAction hq hπmax f hf0 hf1 hf a)
    (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem)
    (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
      (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
      (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
        𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
      (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
      (map_residue_map_twist (integerRingEquiv σ) hf) n (Φ x)
      (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
        f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
        (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hx)
      (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x)))
  rw [map_LubinTateAction hq (integerRingEquiv σ) hπmax hπne0 f hf0 hf1 hf a] at hkey
  exact congrArg (fun w : adjoinIntegers K (Φ x) =>
    ((w : IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure)) hkey

/-! ## 8. ★★★★★払い出し B —— (ii) の有限段 `ρ_{f^σ,n}(Φ τ Φ^{-1}) = σ(ρ_{f,n}(τ))` -/

def reciprocityMap_semilinear_conj.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**申し送りの (ii) の有限段**:

  `ρ_{f^σ, n}(Φ τ Φ^{-1}) = σ(ρ_{f, n}(τ))`。

原文 (Yoshida08 p.9):
> do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

★★**原典と同じ主張ではない**(逸脱として記録する)。原典 Corollary 4.9 の証明の 1 行

> If σ(α) = [xπ_j](α) for σ ∈ W(K^m_f/K), then σ([θ](α)) = [θ]^{(j)}[xπ_j](α) = [xπ′_j][θ](α)
> by Lemma 4.5, hence ρ_{f,m} = ρ_{f′,m}.

は `[θ]`(`𝒪`-**線型**な輸送、`f` と `f′` は**同じ** `K` の上)による共役であるのに対し、
ここでは `Φ`(体の**半線型**自己同型)による共役を扱い、係数側にも `σ` が現れる
(原典で `[θ]^{(j)}` の肩に載っているひねりが、こちらでは `σ` そのものになる)。
★役割(「`ρ` が `f` の取り替えの下でどう動くか」の段)が同じなので `.src` を
Corollary 4.9 に付けた。

★★**`u` の代表元を明示的に取る形**にした(`Units.map` を statement に書かない):
`u` は `ρ_{f,n}(τ)` の代表元、`v` はその `σ` による像の代表元である。
★これで #205(`⇑↑Φ` と `⇑Φ` の食い違い)を statement の段階で避けられる。

★★**証明に使うのは 3 つだけ**:
`reciprocityMap_spec`(在庫、`ρ` の定義性質)、
`existsUnique_unitActionQuotient_eq_algEquiv`(在庫、一意性)、
`semilinear_lubinTateActionAtTorsionPoint`(本ファイル §7)。
★申し送りの **2**(`reciprocityUnits` の spec)と **3**(生成元列の付け替え)は
**有限段では要らなかった** —— 有限段の `ρ_{f,n}` は「点 `x` の上での記述」そのもので、
`x` は任意の `ψ_n` の根でよいからである。 -/
theorem reciprocityMap_semilinear_conj (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a))
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (τ : K.closure ≃ₐ[K.carrier] K.closure) (u v : (𝒪[K.carrier])ˣ)
    (hu : (QuotientGroup.mk u : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n)
      = reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ)
    (hv : (v : 𝒪[K.carrier]) = integerRingEquiv σ (u : 𝒪[K.carrier])) :
    reciprocityMap K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
        (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
        (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (integerRingEquiv σ) hf) n hn (Φ x)
        (map_mem_iteratedLubinTatePsiTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
          f hf0 hf1 hf n hn (Φ : K.closure →+* K.closure)
          (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxψ)
        (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
          f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
          (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxn)
        (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x))
        (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ)
      = QuotientGroup.mk v := by
  have hτx : ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
        (u : 𝒪[K.carrier]) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        : K.closure) = τ x := by
    have h := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ
    rw [← hu] at h
    exact h
  refine (existsUnique_unitActionQuotient_eq_algEquiv K hq
      (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
      (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
      (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
        𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
      (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
      (map_residue_map_twist (integerRingEquiv σ) hf) n hn (Φ x)
      (map_mem_iteratedLubinTatePsiTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
        f hf0 hf1 hf n hn (Φ : K.closure →+* K.closure)
        (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxψ)
      (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
        f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
        (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxn)
      (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x))
      (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ)).unique
    (reciprocityMap_spec K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
      (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
      (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
        𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
      (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
      (map_residue_map_twist (integerRingEquiv σ) hf) n hn (Φ x)
      (map_mem_iteratedLubinTatePsiTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
        f hf0 hf1 hf n hn (Φ : K.closure →+* K.closure)
        (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxψ)
      (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
        f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
        (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxn)
      (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x))
      (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ)) ?_
  show ((unitActionQuotientLift K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
        (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
        (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (integerRingEquiv σ) hf) n (Φ x)
        (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
          f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
          (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxn)
        (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x)) (QuotientGroup.mk v) :
      IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure)
      = conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ (Φ x)
  calc ((unitActionQuotientLift K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
        (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
        (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (integerRingEquiv σ) hf) n (Φ x)
        (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
          f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
          (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxn)
        (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x)) (QuotientGroup.mk v) :
      IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure)
      = ((lubinTateActionAtTorsionPoint K hq
            (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
            (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
            (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
              𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
            (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
            (map_residue_map_twist (integerRingEquiv σ) hf) n (Φ x)
            (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
              f hf0 hf1 hf n (Φ : K.closure →+* K.closure)
              (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ) hxn)
            (IntermediateField.mem_adjoin_simple_self K.carrier (Φ x)) (v : 𝒪[K.carrier]) :
          IntermediateField.adjoin K.carrier ({Φ x} : Set K.closure)) : K.closure) := rfl
    _ = Φ ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
            (u : 𝒪[K.carrier]) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
          : K.closure) := by
        rw [hv]
        exact (semilinear_lubinTateActionAtTorsionPoint K hq σ Φ hΦ hπmax hπne0 f hf0 hf1 hf
          n x hxn hmem (u : 𝒪[K.carrier])).symm
    _ = Φ (τ x) := by rw [hτx]
    _ = conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ (Φ x) := by
        rw [conjSemilinearAlgEquiv_apply, Φ.symm_apply_apply]

/-! ## 9. ★★★`Γ_F` の層への橋

★★直前の波の `absGalConjCME_apply`(`Ψ_g` は `Φ_g` による共役**そのもの**、`rfl`)を
本ファイルの抽象核 C の言葉に翻訳する。★これで §8 の (ii) が
`Γ_F` の `g` に対してそのまま代入できる形になる。 -/

/-- ★★★**`Ψ_g = conjSemilinearAlgEquiv Φ_g`** —— `Γ_F` の共役作用
(`absGalConjCME`)は、本ファイルの抽象核 C(半線型な共役)そのものである。

★両辺とも `y ↦ Φ_g (τ (Φ_g^{-1} y))` に `rfl` で潰れる(#59 の定型 (d):
`closureEquivFixedField` が `RingEquiv` に潰してあるので中間体の層が無い)。

★★**これが「`Γ_F` の層」と §8 の (ii) を繋ぐ橋である。**次の節点は
`reciprocityMap_semilinear_conj` に
`Φ := fixedFieldClosureAut F S hopen g`、
`σ := (fixedFieldAut S g).restrictScalars ℚ_[p]`、
`hΦ := fixedFieldClosureAut_semilinear F S hopen g` を代入し、
この補題で `Ψ_g` に書き換えればよい。 -/
theorem absGalConjCME_eq_conjSemilinearAlgEquiv (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (τ : (fixedFieldLocalField F S hopen).absGal) :
    absGalConjCME F S hopen g τ
      = conjSemilinearAlgEquiv (k := (fixedFieldLocalField F S hopen).carrier)
          (fixedFieldClosureAut F S hopen g)
          ((fixedFieldAut S g).restrictScalars ℚ_[p]).toRingEquiv
          (fixedFieldClosureAut_semilinear F S hopen g) τ :=
  AlgEquiv.ext (fun _ => rfl)

/-! ## 10. ★★★★★副産物 —— **点をまたぐ** Galois 同変性

★★これは `σ := id` の場合(`Φ` が `K.carrier`-線型)で、
★**木が 20 波以上にわたって避けてきた形そのもの**である。

`LubinTateActionEquivariance.lean` は
> `adjoinIntegers K x` と `adjoinIntegers K (σx)` という異なる 2 つの環を
> 橋渡しする必要(cross-point instance bridging)を、`σ(x)` が `K⟮x⟯` に
> 留まるという事実を使って**回避**する

と書き、`lubinTateActionAtAlgEquivPoint`(`x` 自身の座標系で計算した値)を
右辺に置く形でしか同変性を持っていなかった。
★★**本ファイルの §5 を使うと、回避せずに直接書ける。** -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★**点をまたぐ Galois 同変性**:`τ(a · x) = a · (τ x)`。

★右辺は **`adjoinIntegers K (τ x)` の中で計算した本物の Lubin-Tate 作用**である
(`lubinTateActionAtAlgEquivPoint` のような「`x` の座標系での代用品」ではない)。

★★これが「cross-point instance bridging」の**正面突破**である。要ったのは
`adjoinIntegersSemilinear`(§5)だけで、`τ x ∈ K⟮x⟯` は 1 度も使っていない
——実際、この形は `τ x ∉ K⟮x⟯` でも成り立つ(そのときも `K⟮τ x⟯` の中で計算する)。

★証明は §7 と同じで、`σ := AlgEquiv.refl` を入れて係数のひねりを消すだけ
(`PowerSeries.map (RingHom.id) = id`)。 -/
theorem algEquiv_lubinTateActionAtTorsionPoint_cross (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (τ : K.closure ≃ₐ[K.carrier] K.closure) (a : 𝒪[K.carrier]) :
    τ ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)
      = ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n (τ x)
            (algEquiv_mem_iteratedLubinTateTorsionPoints_of_mem K hq hπmax hπne0 f hf0 hf1 hf
              n τ x hx)
            (IntermediateField.mem_adjoin_simple_self K.carrier (τ x)) a :
          IntermediateField.adjoin K.carrier ({τ x} : Set K.closure)) : K.closure) := by
  have hΦ : ∀ c : K.carrier, τ.toRingEquiv (algebraMap K.carrier K.closure c)
      = algebraMap K.carrier K.closure
        ((AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier) c) := fun c => τ.commutes c
  have hid : ∀ c : 𝒪[K.carrier],
      integerRingEquiv (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier) c = c :=
    fun c => Subtype.ext rfl
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  haveI := completeSpace_adjoinIntegers K (τ x)
  haveI := isLinearTopology_adjoinIntegers K (τ x)
  haveI := continuousSMul_adjoinIntegers K (τ x)
  have hkey := aeval_powerSeries_comm_twist (RingHom.id 𝒪[K.carrier])
    ((adjoinIntegersSemilinear K AlgEquiv.refl τ.toRingEquiv hΦ x :
        (adjoinIntegers K x : Type _) ≃+* (adjoinIntegers K (τ x) : Type _)) :
      (adjoinIntegers K x : Type _) →+* (adjoinIntegers K (τ x) : Type _))
    (continuous_adjoinIntegersSemilinear K AlgEquiv.refl τ.toRingEquiv hΦ x)
    (fun c => by
      show adjoinIntegersSemilinear K AlgEquiv.refl τ.toRingEquiv hΦ x
          (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) c)
        = algebraMap 𝒪[K.carrier] (adjoinIntegers K (τ x)) c
      rw [semilinear_adjoinIntegersSemilinear, hid]
      rfl)
    (LubinTateAction hq hπmax f hf0 hf1 hf a)
    (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem)
    (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n (τ x)
      (algEquiv_mem_iteratedLubinTateTorsionPoints_of_mem K hq hπmax hπne0 f hf0 hf1 hf n τ x hx)
      (IntermediateField.mem_adjoin_simple_self K.carrier (τ x)))
  simp only [PowerSeries.map_id, id_eq] at hkey
  exact congrArg (fun w : adjoinIntegers K (τ x) =>
    ((w : IntermediateField.adjoin K.carrier ({τ x} : Set K.closure)) : K.closure)) hkey

end ABC3.Found.PGC
