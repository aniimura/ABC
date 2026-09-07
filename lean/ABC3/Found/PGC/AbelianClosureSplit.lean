import ABC3.Found.PGC.LubinTateZhat
import ABC3.Found.PGC.ArithFrobeniusTopGen
import ABC3.Found.PGC.AbelianFrobeniusSplit
import ABC3.Found.PGC.TopAbelianization

/-!
# Frobenius 持ち上げ `σ` と `K^ab = K^ur·E_σ` —— [Yoshida08] Theorem 6.15 の B2・B3

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

原文の Proof 段落(`0_Source` の `.txt` 1180-1200 行。◆抽出器の字形の潰れは直していない):

> Proof. Take a σ ∈W(KLT/K) with v(σ) = n > 0, and let L = Kn. Extend σ arbitrarily
> to σ ∈W(Kab/K), and let Eσ ⊂Kab be its fixed field. Then Eσ ∩Kur = L and Eσ/L is
> totally ramified Galois. Now Gal(Kab/Eσ) ∼= bZ with σ 7→1 by the definition of Eσ. On the
> other hand, Gal(KurEσ/Eσ) ∼= Gal(Kur/L) ∼= bZ by σ 7→1, as σ|Kur = FrobL. Therefore
> Gal(Kab/Eσ) ∼= Gal(KurEσ/Eσ), i.e. Kab = KurEσ.

## ★本ノードの範囲

★**本ノードは Thm 6.15 の B2・B3 のみ**である。決定 D29 により、原文の
「Take a σ with v(σ) = n > 0」を **`n = 1` に固定**する
(原文は「任意の σ について」ではなく「ある σ を取れ」と言っているので、
これは逸脱ではなく選択の固定である)。`n = 1` では `L = K` になり、
相対 Lubin-Tate が絶対 Lubin-Tate に潰れて木の機構の字面と一致する。

* **B2**: `σ ∈ Γ_K` で `σ|_{K_π} = id` かつ `σ|_{K^ur} = 算術 Frobenius`
  (`Gal(K^ur/K) ≅ Ẑ` の位相的生成元)なるものの存在。
* **B3**: その `σ` について `Ω = K^ur·(Ω ⊓ E_σ)`(`E_σ = K̄^σ`)。
  ★`Ω := K^ab` と取れば原文の `K^ab = K^ur E_σ` そのものである。

B1・B4・B5・B6(`K^LT ≤ K^ab`、Prop 6.14 の n=1 形、`E_σ ⊆ K_π`、組み立て)は
本ノードの担当ではない。

## ★`K^ab` について(測定の記録)

`node tools/decl-index.mjs` で作った `.cache/decl-index.txt`(24,942 宣言)を
`grep -n "abelianClosure\|maximalAbelian"` で引いたところ **0 件**、
`node tools/absent-recheck.mjs --try 'abelian.*[Cc]losure|maximalAbelian|K\^ab'` も
ABC3 側 0 件であった。すなわち**この木には `K^ab` という名前がまだ無い**。

そこで B3 は

```
U ⊔ (Ω ⊓ fixedField ⟨σ⟩) = Ω        (U = K^ur ≤ Ω、Ω は K 上正規)
```

という **`Ω` を動かせる形**で証明し、次の 2 つを系として出してある:

* `Ω := ⊤`(= `K̄`)—— `K̄ = K^ur·K̄^σ`
* `Ω := IntermediateField.fixedField ((commutator Γ_K).topologicalClosure)`
  —— ★★**これが `K^ab` であり、原文の `K^ab = K^ur E_σ` そのもの**

★★**新しい `def` は置いていない。** `K^ab` は上の式で**書き下して**ある。
専用ノードが `abelianClosure` を立てたときに**同名衝突**
(`lean-idioms.md` #146:同名が 2 つあると曖昧参照で落ちる)を起こさないためである。
立ったら定義が `rfl` なので、そのまま本節の主張に乗り移れる。

★`K^ur ≤ K^ab` は `commutator_le_fixingSubgroup_unramifiedClosure`
(= `K^ur/K` はアーベル)として本ファイルで証明した。材料は
`isCyclic_gal_of_isUnramifiedAdjoin`(`Found/PGC/AbelianFrobeniusSplit.lean`)と
`mem_unramifiedClosure_iff`(`Found/PGC/UnramifiedExtension.lean`)。

## ★退化の自己検査

* ★★**`E_σ` は省けない。** Corollary 6.13(ii) は `K′K′′/K` が完全分岐であることを
  仮定するが、「2 つの完全分岐アーベル拡大の合成」で代用すると**壊れる**:
  `K = ℚ_p` で `ℚ_p(√p)` と `ℚ_p(√(up))`(`u` は非平方単数)の合成は
  `√u` を含み、`ℚ_p(√u)/ℚ_p` は**不分岐**である。完全分岐性は合成で保たれない。
  ゆえに完全分岐性は `E_σ`(完全分岐)への包含から供給するほかなく、
  B2(`σ` の構成)と B3(`K^ab = K^ur E_σ`)を飛ばすことはできない。
* ★**`σ|_{K_π} = id` と `σ|_{K^ur} = Frobenius` の両方が要る。**
  前者を落とすと `E_σ` が小さすぎて `E_σ ⊆ K_π`(B5)が出ず、
  後者を落とすと `E_σ` が大きすぎて `K^ur·E_σ = K^ab` が破れる
  (例:`σ = 1` なら `E_σ = K^ab` で `K^ur·E_σ = K^ab` は自明に成り立つが、
  `E_σ ⊆ K_π` が偽になる)。
* ★**示すのは `K^ab = K^ur E_σ` であって `K^ab = K^LT` ではない。**
  後者は B5(`E_σ ⊆ K_π`)・B6 を経由して初めて出る。
* 抽象核 `mem_of_mem_closure_zpowers` の `hkey`(`σ^k ∈ C → e ∣ k`)は落とせない。
  落とすと `C = ⊤` が反例になる(`τ` が `N` に入る保証が消える)。
* 抽象核 `sup_inf_fixedField_zpowers_eq` の `hUΩ : U ≤ Ω` は落とせない。
  落とすと左辺が `Ω` を超える。

## ☆逸脱の記録

* D29 により `n = 1` に固定した(上記)。原文の `L = K_n` は `L = K` になる。
* 原文は `E_σ ⊂ K^ab` を「`K^ab` の中の `σ` の固定体」として導入するが、
  本ファイルは `K̄` の中の固定体 `K̄^σ` を使い、`Ω ⊓ K̄^σ` として原文の `E_σ` を
  切り出す。`Ω = K^ab` では両者は一致する(`K^ab` の中の `σ|_{K^ab}` の固定体
  = `K^ab ∩ K̄^σ`)。★これは `K^ab` がまだ木に無いための読み替えであり、
  後続の消費(B5・B6)には影響しない。
* 原文は `Gal(K^ab/E_σ) ≅ Ẑ ≅ Gal(K^ur E_σ/E_σ)` という **`Ẑ` の Hopf 性**
  (副有限群の全射自己準同型は単射)を経由するが、本ファイルは
  **`Ẑ` を一切使わない**初等的な道を取った:
  有限次 Galois 中間体 `M` ごとに `e := ord(σ|_M)` を取り、
  不分岐塔の段 `K_e` を噛ませて `e ∣ k` を絞る。
  ★★**原典より短い道**である(`Ẑ` の Hopf 性の形式化が丸ごと不要になった)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## 0. 抽象核 A —— 位相群だけの言葉

★分岐・付値・Galois の語彙が 1 つも出てこない。 -/

/-- ★★**抽象核 A**:`τ` が `⟨σ⟩` の閉包に属し、`C` が「`σ^k ∈ C ⟹ e ∣ k`」を
検出する開部分群で `τ ∈ C` なら、`σ^e ∈ N` なる任意の開部分群 `N` に `τ` が入る。

`τ` の開近傍 `τ·(N ⊓ C)` が `⟨σ⟩` と交わることから `σ^k = τ n` と書き、
`σ^k ∈ C`(`τ`・`n` がともに `C` に居る)⟹ `e ∣ k` ⟹ `σ^k ∈ N`、
最後に `τ = σ^k n⁻¹ ∈ N`。

★`hkey` を落とすと偽(`C = ⊤` が反例)。 -/
theorem mem_of_mem_closure_zpowers {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] (σ τ : G) (N C : Subgroup G) (e : ℕ)
    (hτ : τ ∈ closure ((Subgroup.zpowers σ : Subgroup G) : Set G))
    (hNopen : IsOpen (N : Set G)) (hCopen : IsOpen (C : Set G))
    (hτC : τ ∈ C) (he : σ ^ e ∈ N)
    (hkey : ∀ k : ℤ, σ ^ k ∈ C → (e : ℤ) ∣ k) :
    τ ∈ N := by
  set U : Set G := (fun g => τ * g) '' ((N ⊓ C : Subgroup G) : Set G) with hU
  have hopen : IsOpen U := isOpenMap_mul_left τ _ (hNopen.inter hCopen)
  have hτU : τ ∈ U := ⟨1, ⟨N.one_mem, C.one_mem⟩, by simp⟩
  obtain ⟨y, hyU, hyz⟩ := mem_closure_iff.mp hτ U hopen hτU
  obtain ⟨n, ⟨hnN, hnC⟩, hn⟩ := hyU
  obtain ⟨k, hk⟩ := hyz
  have hn' : τ * n = y := hn
  have hk' : σ ^ k = y := hk
  have hkC : σ ^ k ∈ C := by
    rw [hk', ← hn']
    exact C.mul_mem hτC hnC
  obtain ⟨m, hm⟩ := hkey k hkC
  have hkN : σ ^ k ∈ N := by
    rw [hm, zpow_mul, zpow_natCast]
    exact N.zpow_mem he m
  have hτeq : τ = σ ^ k * n⁻¹ := by rw [hk', ← hn']; group
  rw [hτeq]
  exact N.mul_mem hkN (N.inv_mem hnN)

/-- ★**抽象核 A'**:コンパクト Hausdorff 位相群では、閉部分群と**正規**閉部分群の
`⊔` は閉。`↑(A ⊔ J) = ↑A * ↑J`(`Subgroup.mul_normal`)と、コンパクト集合の積が
コンパクトであることから。 -/
theorem isClosed_sup_of_normal {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [T2Space G] (A J : Subgroup G) [J.Normal]
    (hA : IsCompact (A : Set G)) (hJ : IsCompact (J : Set G)) :
    IsClosed ((A ⊔ J : Subgroup G) : Set G) := by
  rw [Subgroup.mul_normal]
  exact (hA.mul hJ).isClosed

/-- `J` が正規なら `A ⊔ J` の元は `a * j` と書ける。 -/
theorem exists_mul_of_mem_sup_normal {G : Type*} [Group G] (A J : Subgroup G) [J.Normal]
    {τ : G} (h : τ ∈ A ⊔ J) : ∃ a ∈ A, ∃ j ∈ J, a * j = τ := by
  have hmem : τ ∈ ((A ⊔ J : Subgroup G) : Set G) := h
  rw [Subgroup.mul_normal] at hmem
  obtain ⟨a, ha, j, hj, haj⟩ := hmem
  exact ⟨a, ha, j, hj, haj⟩

/-! ## 1. 抽象核 B —— 無限次 Galois 対応だけの言葉

★ここも分岐・付値・Lubin-Tate の語彙は出てこない。 -/

section AbstractGalois

variable {k E : Type*} [Field k] [Field E] [Algebra k E]

/-- ★★**抽象核 B**:`Gal(E/A) ∩ Gal(E/B) = 1` なら `A ⊔ B = ⊤`。

`(A ⊔ B).fixingSubgroup ≤ A.fixingSubgroup ⊓ B.fixingSubgroup = ⊥` から
`A ⊔ B = fixedField ⊥ = ⊤`(`InfiniteGalois.fixedField_fixingSubgroup`)。 -/
theorem sup_eq_top_of_forall_eq_one [IsGalois k E] (A B : IntermediateField k E)
    (h : ∀ σ : E ≃ₐ[k] E, σ ∈ A.fixingSubgroup → σ ∈ B.fixingSubgroup → σ = 1) :
    A ⊔ B = ⊤ := by
  have hbot : (A ⊔ B).fixingSubgroup = ⊥ := by
    refine le_antisymm (fun σ hσ => ?_) bot_le
    exact (Subgroup.mem_bot).mpr
      (h σ (IntermediateField.fixingSubgroup_le (le_sup_left : A ≤ A ⊔ B) hσ)
        (IntermediateField.fixingSubgroup_le (le_sup_right : B ≤ A ⊔ B) hσ))
  have hff := InfiniteGalois.fixedField_fixingSubgroup (k := k) (K := E) (A ⊔ B)
  rw [hbot, IntermediateField.fixedField_bot] at hff
  exact hff.symm

/-- ★**抽象核 C**:`fixedField H` を各点固定する元は `H` の**位相的閉包**に入る。

`InfiniteGalois.fixingSubgroup_fixedField` を閉部分群 `H̄` に当て、
`fixedField H̄ ≤ fixedField H`(`H ≤ H̄` の反変)で引き戻す。 -/
theorem fixingSubgroup_fixedField_le_topologicalClosure [IsGalois k E]
    (H : Subgroup (E ≃ₐ[k] E)) :
    (IntermediateField.fixedField H).fixingSubgroup ≤ H.topologicalClosure := by
  have hcl : (IntermediateField.fixedField (H.topologicalClosure)).fixingSubgroup
      = H.topologicalClosure :=
    InfiniteGalois.fixingSubgroup_fixedField
      ⟨H.topologicalClosure, H.isClosed_topologicalClosure⟩
  rw [← hcl]
  refine IntermediateField.fixingSubgroup_le ?_
  intro x hx
  rw [IntermediateField.mem_fixedField_iff] at hx ⊢
  exact fun g hg => hx g (H.le_topologicalClosure hg)

/-- ★★**抽象核 D**(B2 の核):`A ⊓ B = k` のとき、`Gal(B/k)` の任意の元 `b` は
「**`A` を各点固定する**」`σ ∈ Gal(E/k)` に持ち上がる。

`restrictPairHom_surjective`(`Found/PGC/AbelianDecomposition.lean`)に
`(1, b)` を渡すだけ。★`A ⊓ B = ⊥` を落とすと偽
(`A = B ≠ ⊥` なら `b ≠ 1` を `A` を固定する元に持ち上げられない)。 -/
theorem exists_mem_fixingSubgroup_restrictNormalHom_eq [IsGalois k E]
    (A B : IntermediateField k E) [Normal k A] [Normal k B] (h : A ⊓ B = ⊥)
    (b : (B : Type _) ≃ₐ[k] (B : Type _)) :
    ∃ σ : E ≃ₐ[k] E, σ ∈ A.fixingSubgroup ∧
      AlgEquiv.restrictNormalHom (F := k) (K₁ := E) ((B : IntermediateField k E) : Type _) σ = b := by
  obtain ⟨σ, hσ⟩ := restrictPairHom_surjective A B h (1, b)
  refine ⟨σ, ?_, ?_⟩
  · exact (restrictNormalHom_eq_one_iff' A σ).mp (congrArg Prod.fst hσ)
  · exact congrArg Prod.snd hσ

/-- ★★★**抽象核 E**(B3 の核):`U ≤ Ω`(`Ω` は `k` 上正規)で、
「`U` を各点固定し、かつ `⟨σ⟩` の閉包に属する元は `1` に限る」なら

```
U ⊔ (Ω ⊓ E_σ) = Ω        (E_σ := fixedField ⟨σ⟩)
```

★★**原典の `Ẑ` の Hopf 性を使わない。** 使うのは
`fixedField (J ⊔ ⟨σ⟩) = Ω ⊓ E_σ`(`fixedField_sup`)と、
`J` が正規閉なので `⟨σ⟩‾ ⊔ J` が閉であること(抽象核 A')だけ。

★分岐・付値・Lubin-Tate の語彙は 1 つも出てこない。 -/
theorem sup_inf_fixedField_zpowers_eq [IsGalois k E] (U Ω : IntermediateField k E)
    [Normal k Ω] (σ : E ≃ₐ[k] E) (hUΩ : U ≤ Ω)
    (hkey : ∀ τ : E ≃ₐ[k] E, τ ∈ U.fixingSubgroup →
      τ ∈ (Subgroup.zpowers σ).topologicalClosure → τ = 1) :
    U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) = Ω := by
  haveI : (Ω.fixingSubgroup).Normal := normal_fixingSubgroup Ω
  have hXΩ : U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) ≤ Ω :=
    sup_le hUΩ inf_le_left
  have hfix : (U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))).fixingSubgroup
      = Ω.fixingSubgroup := by
    refine le_antisymm ?_ (IntermediateField.fixingSubgroup_le hXΩ)
    intro τ hτ
    have hτU : τ ∈ U.fixingSubgroup :=
      IntermediateField.fixingSubgroup_le
        (le_sup_left : U ≤ U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))) hτ
    have hτI : τ ∈ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)).fixingSubgroup :=
      IntermediateField.fixingSubgroup_le
        (le_sup_right : (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) ≤
          U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))) hτ
    have heq : Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)
        = IntermediateField.fixedField (Ω.fixingSubgroup ⊔ Subgroup.zpowers σ) := by
      rw [fixedField_sup, InfiniteGalois.fixedField_fixingSubgroup]
    rw [heq] at hτI
    have hτcl : τ ∈ (Ω.fixingSubgroup ⊔ Subgroup.zpowers σ).topologicalClosure :=
      fixingSubgroup_fixedField_le_topologicalClosure _ hτI
    have hle : (Ω.fixingSubgroup ⊔ Subgroup.zpowers σ).topologicalClosure
        ≤ (Subgroup.zpowers σ).topologicalClosure ⊔ Ω.fixingSubgroup := by
      refine Subgroup.topologicalClosure_minimal _ (sup_le ?_ ?_) ?_
      · exact le_sup_right
      · exact le_trans (Subgroup.le_topologicalClosure _) le_sup_left
      · exact isClosed_sup_of_normal _ _
          (Subgroup.isClosed_topologicalClosure _).isCompact
          (InfiniteGalois.fixingSubgroup_isClosed Ω).isCompact
    obtain ⟨z, hz, j, hj, hzj⟩ := exists_mul_of_mem_sup_normal _ _ (hle hτcl)
    have hjU : j ∈ U.fixingSubgroup := IntermediateField.fixingSubgroup_le hUΩ hj
    have hzU : z ∈ U.fixingSubgroup := by
      have : z = τ * j⁻¹ := by rw [← hzj]; group
      rw [this]
      exact Subgroup.mul_mem _ hτU (Subgroup.inv_mem _ hjU)
    have hz1 : z = 1 := hkey z hzU hz
    have : τ = j := by rw [← hzj, hz1, one_mul]
    rw [this]
    exact hj
  calc U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))
      = IntermediateField.fixedField
          (U ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))).fixingSubgroup :=
        (InfiniteGalois.fixedField_fixingSubgroup _).symm
    _ = IntermediateField.fixedField Ω.fixingSubgroup := by rw [hfix]
    _ = Ω := InfiniteGalois.fixedField_fixingSubgroup Ω

end AbstractGalois

/-! ## 2. 具体層 B2 —— `K_π` を各点固定する Frobenius 持ち上げ -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★**B2(一般形)**:`F ⊓ K^ur = K` なる `K` 上正規な中間体 `F` について、
`F` を各点固定し、`K^ur` 上は**算術 Frobenius**(`Gal(K^ur/K) ≅ Ẑ` の位相的生成元)
になる `σ ∈ Γ_K` が在る。

抽象核 D に `A := F`、`B := K^ur`、`b := arithFrobenius K` を代入するだけ。

★退化:`F ⊓ K^ur = ⊥` は落とせない。例えば `F = K_2`(2 次不分岐)とすると、
`F` を固定する `σ` の `Gal(K_2/K)` への像は `1` で、Frobenius になりえない。 -/
theorem exists_arithFrobenius_lift_fixing (K : PAdicLocalField p)
    (F : IntermediateField K.carrier K.closure) (hFn : Normal K.carrier F)
    (hF : F ⊓ unramifiedClosure K = ⊥) :
    ∃ σ : K.absGal, σ ∈ F.fixingSubgroup ∧
      AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        (↥(unramifiedClosure K)) σ = arithFrobenius K := by
  haveI := isGalois_closure K
  haveI := hFn
  haveI := normal_unramifiedClosure K
  exact exists_mem_fixingSubgroup_restrictNormalHom_eq F (unramifiedClosure K) hF
    (arithFrobenius K)

/-- ★★★**B2**:`K̄` の中に `K` 上正規な `K_π`(`Gal(K_π/K) ≅ 𝒪_K^×`)と
`σ ∈ Γ_K` が在って

* `K_π ⊓ K^ur = K`
* `σ|_{K_π} = id`
* `σ|_{K^ur} = 算術 Frobenius`

★決定 D29 により `n = 1`(原文の `L = K_n` は `L = K`)。
★`K_π` は `∃` の内側に閉じ込めてある。 -/
theorem exists_lubinTateClosure_arithFrobenius_lift (K : PAdicLocalField p) :
    ∃ (E : IntermediateField K.carrier K.closure) (σ : K.absGal),
      Normal K.carrier E ∧ E ⊓ unramifiedClosure K = ⊥ ∧
      Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) ∧
      σ ∈ E.fixingSubgroup ∧
      AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        (↥(unramifiedClosure K)) σ = arithFrobenius K := by
  obtain ⟨E, hEn, hEinf, hEgal, -⟩ := exists_lubinTateUnramified_decomposition K
  obtain ⟨σ, hσF, hσur⟩ := exists_arithFrobenius_lift_fixing K E hEn hEinf
  exact ⟨E, σ, hEn, hEinf, hEgal, hσF, hσur⟩

/-! ## 3. 具体層 B3 —— `⟨σ⟩` の閉包は惰性群と自明にしか交わらない -/

/-- ★★★**B3 の心臓部**:`σ|_{K^ur} = 算術 Frobenius` のとき、
`K^ur` を各点固定し、かつ `⟨σ⟩` の閉包に属する `τ ∈ Γ_K` は `1` に限る。
すなわち `Gal(K̄/K^ur) ⊓ ⟨σ⟩‾ = 1`。

## ★原典より短い道

原文は `Gal(K^ab/E_σ) ≅ Ẑ ≅ Gal(K^ur E_σ/E_σ)` という **`Ẑ` の Hopf 性**を経由するが、
ここでは `Ẑ` を一切使わない:

`x ∈ K̄` を取り、`x` を含む有限次 Galois 中間体 `M`
(`FiniteGaloisIntermediateField.adjoin`)と `e := ord(σ|_M) ≠ 0` を置く。
不分岐塔の段 `K_e` の固定部分群の引き戻しを `C` とすると

* `τ ∈ C`(`τ` は `K^ur` を固定するので `τ|_{K^ur} = 1`)
* `σ^k ∈ C ⟹ e ∣ k`(`zpow_arithFrobenius_mem_fixingSubgroup_iff`)
* `σ^e ∈ Gal(K̄/M)`(`e` は `σ|_M` の位数)

なので抽象核 A が `τ ∈ Gal(K̄/M)`、したがって `τ x = x` を与える。

★★`hσ` は落とせない。`σ = 1` とすると `⟨σ⟩‾ = 1` で主張は自明になるが、
`σ|_{K^ur}` が Frobenius の**べき**(位数有限の像)だと `e` を検出できず偽になる
(例:`σ ∈ Gal(K̄/K^ur)` が無限位数なら `⟨σ⟩‾ ⊆ Gal(K̄/K^ur)` で反例)。 -/
theorem eq_one_of_fixes_unramifiedClosure_of_mem_closure_zpowers (K : PAdicLocalField p)
    {σ : K.absGal}
    (hσ : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    {τ : K.absGal} (hτur : τ ∈ (unramifiedClosure K).fixingSubgroup)
    (hτ : τ ∈ (Subgroup.zpowers σ).topologicalClosure) :
    τ = 1 := by
  haveI := isGalois_closure K
  haveI := normal_unramifiedClosure K
  have hres1 : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) τ = 1 :=
    (restrictNormalHom_eq_one_iff' (unramifiedClosure K) τ).mpr hτur
  have key : ∀ x : K.closure, τ x = x := by
    intro x
    let M : IntermediateField K.carrier K.closure :=
      (FiniteGaloisIntermediateField.adjoin K.carrier ({x} : Set K.closure)).toIntermediateField
    have hxM : x ∈ M :=
      FiniteGaloisIntermediateField.subset_adjoin K.carrier ({x} : Set K.closure) rfl
    haveI : FiniteDimensional K.carrier M := inferInstance
    haveI : Normal K.carrier M := inferInstance
    haveI : Finite ((M : Type _) ≃ₐ[K.carrier] (M : Type _)) := inferInstance
    obtain ⟨e, he0, heM⟩ : ∃ e : ℕ, e ≠ 0 ∧ σ ^ e ∈ M.fixingSubgroup := by
      refine ⟨orderOf (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        ((M : IntermediateField K.carrier K.closure) : Type _) σ), (orderOf_pos _).ne', ?_⟩
      rw [← IntermediateField.restrictNormalHom_ker M, MonoidHom.mem_ker, map_pow]
      exact pow_orderOf_eq_one _
    refine (IntermediateField.mem_fixingSubgroup_iff M τ).mp ?_ x hxM
    refine mem_of_mem_closure_zpowers σ τ M.fixingSubgroup
      (((unramLevel K e).fixingSubgroup).comap
        (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
          (↥(unramifiedClosure K)))) e hτ
      (IntermediateField.fixingSubgroup_isOpen M) ?_ ?_ heM ?_
    · exact (isOpen_fixingSubgroup_unramLevel K e).preimage
        (InfiniteGalois.restrictNormalHom_continuous (unramifiedClosure K))
    · refine Subgroup.mem_comap.mpr ?_
      rw [hres1]
      exact one_mem _
    · intro k hk
      have hk' := Subgroup.mem_comap.mp hk
      rw [map_zpow, hσ] at hk'
      exact (zpow_arithFrobenius_mem_fixingSubgroup_iff K he0 k).mp hk'
  refine AlgEquiv.ext fun x => ?_
  rw [key x, AlgEquiv.one_apply]

/-- ★★★**B3(`Ω` 版)**:`σ|_{K^ur} = 算術 Frobenius` のとき、`K^ur` を含む
`K` 上正規な任意の中間体 `Ω` について

```
Ω = K^ur · (Ω ⊓ K̄^σ)
```

★★`Ω := K^ab` と取れば、これが原文の `K^ab = K^ur E_σ` である
(`E_σ = K^ab ⊓ K̄^σ`)。★`abelianClosure` はまだこの木に無いので、
`Ω` を仮定として置いてある —— 立った時点で代入するだけで済む。 -/
theorem unramifiedClosure_sup_inf_fixedField_zpowers_eq (K : PAdicLocalField p)
    (Ω : IntermediateField K.carrier K.closure) (hΩn : Normal K.carrier Ω)
    (hΩ : unramifiedClosure K ≤ Ω) {σ : K.absGal}
    (hσ : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K) :
    unramifiedClosure K ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) = Ω := by
  haveI := isGalois_closure K
  haveI := hΩn
  exact sup_inf_fixedField_zpowers_eq (unramifiedClosure K) Ω σ hΩ
    fun τ h1 h2 => eq_one_of_fixes_unramifiedClosure_of_mem_closure_zpowers K hσ h1 h2

/-- ★★**B3(`Ω = K̄` 版)**:`K̄ = K^ur · K̄^σ`。

★`Ω` 版の系ではなく抽象核 B から直接出す(`Normal K.carrier ⊤` を避けるため)。 -/
theorem unramifiedClosure_sup_fixedField_zpowers_eq_top (K : PAdicLocalField p)
    {σ : K.absGal}
    (hσ : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K) :
    unramifiedClosure K ⊔ IntermediateField.fixedField (Subgroup.zpowers σ) = ⊤ := by
  haveI := isGalois_closure K
  refine sup_eq_top_of_forall_eq_one (unramifiedClosure K)
    (IntermediateField.fixedField (Subgroup.zpowers σ)) fun τ h1 h2 => ?_
  exact eq_one_of_fixes_unramifiedClosure_of_mem_closure_zpowers K hσ h1
    (fixingSubgroup_fixedField_le_topologicalClosure _ h2)

/-! ## 4. `K^ab` への具体化 —— 原文の `K^ab = K^ur E_σ` そのもの

★**`def` を新たに置かない。** `K^ab` は
`IntermediateField.fixedField ((commutator Γ_K).topologicalClosure)` と**書き下す**。
専用ノードが `abelianClosure` を立てたときに**同名衝突**(`lean-idioms.md` #146)を
起こさないためである。立ったら `rfl` で本節の主張に乗り移れる。 -/

/-- 不分岐な単項拡大の Galois 群は巡回(`isCyclic_gal_of_isUnramifiedAdjoin`)、
ゆえに可換。したがって `Γ_K` の交換子は `K(x)` を各点固定する。 -/
theorem commutatorElement_mem_fixingSubgroup_of_isUnramifiedAdjoin (K : PAdicLocalField p)
    {x : K.closure} (hx : IsUnramifiedAdjoin K x) (a b : K.absGal) :
    a * b * a⁻¹ * b⁻¹ ∈
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup := by
  haveI := normal_of_isUnramifiedAdjoin K x hx
  haveI := isCyclic_gal_of_isUnramifiedAdjoin K x hx
  set Kx := IntermediateField.adjoin K.carrier ({x} : Set K.closure) with hKx
  set φ := AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (Kx : Type _) with hφ
  have hker : φ.ker = Kx.fixingSubgroup := IntermediateField.restrictNormalHom_ker Kx
  obtain ⟨g0, hg0⟩ := IsCyclic.exists_generator (α := (Kx ≃ₐ[K.carrier] Kx))
  rw [← hker, MonoidHom.mem_ker]
  simp only [map_mul, map_inv]
  exact commutator_eq_one_of_mul_comm (mul_comm_of_forall_mem_zpowers hg0 (φ a) (φ b))

/-- ★**`⁅Γ_K, Γ_K⁆ ≤ Gal(K̄/K^ur)`**——すなわち `K^ur/K` はアーベル。

`mem_unramifiedClosure_iff`(`K^ur` の元はどれか 1 つの不分岐 `K(x)` に入る)で
有限段に落として上の補題を当てる。★**2 層(`K̄` と `K^ur`)をまたがない**ので
`lean-idioms.md` #59 に触れない。 -/
theorem commutator_le_fixingSubgroup_unramifiedClosure (K : PAdicLocalField p) :
    commutator K.absGal ≤ (unramifiedClosure K).fixingSubgroup := by
  rw [commutator_def]
  refine Subgroup.commutator_le.mpr fun a _ b _ => ?_
  show a * b * a⁻¹ * b⁻¹ ∈ (unramifiedClosure K).fixingSubgroup
  rw [IntermediateField.mem_fixingSubgroup_iff]
  intro z hz
  obtain ⟨x, hx, hzx⟩ := (mem_unramifiedClosure_iff K z).mp hz
  exact (IntermediateField.mem_fixingSubgroup_iff _ _).mp
    (commutatorElement_mem_fixingSubgroup_of_isUnramifiedAdjoin K hx a b) z hzx

/-- ★**`K^ur ≤ K^ab`**。`Gal(K̄/K^ur)` は閉なので位相的閉包でも成り立つ。 -/
theorem unramifiedClosure_le_fixedField_commutator (K : PAdicLocalField p) :
    unramifiedClosure K
      ≤ IntermediateField.fixedField ((commutator K.absGal).topologicalClosure) := by
  haveI := isGalois_closure K
  refine (IntermediateField.le_iff_le _ _).mpr ?_
  exact Subgroup.topologicalClosure_minimal _
    (commutator_le_fixingSubgroup_unramifiedClosure K)
    (InfiniteGalois.fixingSubgroup_isClosed _)

/-- ★**`K^ab/K` は Galois**(したがって `Normal`)。
`(fixedField H̄).fixingSubgroup = H̄` と `⁅Γ,Γ⁆‾` の正規性から。 -/
theorem isGalois_fixedField_commutator (K : PAdicLocalField p) :
    IsGalois K.carrier
      (IntermediateField.fixedField ((commutator K.absGal).topologicalClosure)) := by
  haveI := isGalois_closure K
  have hfix : (IntermediateField.fixedField
        ((commutator K.absGal).topologicalClosure)).fixingSubgroup
      = (commutator K.absGal).topologicalClosure :=
    InfiniteGalois.fixingSubgroup_fixedField
      ⟨(commutator K.absGal).topologicalClosure, Subgroup.isClosed_topologicalClosure _⟩
  rw [← InfiniteGalois.normal_iff_isGalois, hfix]
  infer_instance

/-- ★★★★**B3 の逐語形**:`K^ab = K^ur·E_σ`(`E_σ = K^ab ⊓ K̄^σ`)。

原文:
> Therefore Gal(Kab/Eσ) ∼= Gal(KurEσ/Eσ), i.e. Kab = KurEσ.

★`Ω := K^ab = fixedField (⁅Γ_K,Γ_K⁆‾)` を
`unramifiedClosure_sup_inf_fixedField_zpowers_eq` に代入しただけ。 -/
theorem fixedField_commutator_eq_unramifiedClosure_sup_inf_fixedField_zpowers
    (K : PAdicLocalField p) {σ : K.absGal}
    (hσ : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K) :
    unramifiedClosure K ⊔
        (IntermediateField.fixedField ((commutator K.absGal).topologicalClosure) ⊓
          IntermediateField.fixedField (Subgroup.zpowers σ))
      = IntermediateField.fixedField ((commutator K.absGal).topologicalClosure) := by
  haveI := isGalois_fixedField_commutator K
  exact unramifiedClosure_sup_inf_fixedField_zpowers_eq K _ IsGalois.to_normal
    (unramifiedClosure_le_fixedField_commutator K) hσ

/-! ## 5. B2 + B3 をひとまとめに -/

/-- ★★★★**Theorem 6.15 の B2 + B3**。

`K̄` の中に `K` 上正規な `K_π`(`Gal(K_π/K) ≅ 𝒪_K^×`)と `σ ∈ Γ_K` が在って

* `K_π ⊓ K^ur = K`
* `σ|_{K_π} = id`(B2)
* `σ|_{K^ur} = 算術 Frobenius`(B2)
* `K^ur` を含む `K` 上正規な任意の `Ω` について `Ω = K^ur·(Ω ⊓ K̄^σ)`(B3)
* とくに `Ω := K^ab = fixedField (⁅Γ_K,Γ_K⁆‾)` で
  **`K^ab = K^ur·E_σ`**(原文の逐語形、B3)

★**本ノードは Thm 6.15 の B2・B3 のみ**であり、`K^ab = K^LT` は
B5(`E_σ ⊆ K_π`)・B6 の担当である。 -/
theorem exists_frobeniusLift_abelian_split (K : PAdicLocalField p) :
    ∃ (E : IntermediateField K.carrier K.closure) (σ : K.absGal),
      Normal K.carrier E ∧ E ⊓ unramifiedClosure K = ⊥ ∧
      Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) ∧
      σ ∈ E.fixingSubgroup ∧
      AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        (↥(unramifiedClosure K)) σ = arithFrobenius K ∧
      (∀ Ω : IntermediateField K.carrier K.closure, Normal K.carrier Ω →
        unramifiedClosure K ≤ Ω →
        unramifiedClosure K ⊔ (Ω ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) = Ω) ∧
      unramifiedClosure K ⊔
          (IntermediateField.fixedField ((commutator K.absGal).topologicalClosure) ⊓
            IntermediateField.fixedField (Subgroup.zpowers σ))
        = IntermediateField.fixedField ((commutator K.absGal).topologicalClosure) := by
  obtain ⟨E, σ, hEn, hEinf, hEgal, hσF, hσur⟩ :=
    exists_lubinTateClosure_arithFrobenius_lift K
  exact ⟨E, σ, hEn, hEinf, hEgal, hσF, hσur,
    fun Ω hΩn hΩ => unramifiedClosure_sup_inf_fixedField_zpowers_eq K Ω hΩn hΩ hσur,
    fixedField_commutator_eq_unramifiedClosure_sup_inf_fixedField_zpowers K hσur⟩

def exists_frobeniusLift_abelian_split.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

end ABC3.Found.PGC
