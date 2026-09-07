import ABC3.Found.PGC.CyclotomicRecoverable

/-!
# 経路 Λ12 —— Artin 写像の同変性(壁を `μ_{p^n}(K̄)` の言葉に書き換える)

[pGC] Proposition 1.1 への**経路 Λ**の最後の 1 本。Λ11
(`Found/PGC/CyclotomicRecoverable.lean`)は Proposition 1.1 を

  `CyclotomeConjIsCyclotomic p`(存在量化子を 1 つも含まない 1 文)

に還元した。本ファイルはその 1 文が**何であるか**を確定させる:

  「開正規 `S ⊴ Γ_F` の円分子 `Λ_n(S)` から `μ_{p^n}(K̄)` への
   **`Γ_F`-同変な**群同型が存在する。」

すなわち残った壁は、Artin 写像の同変性
`Art_L(σ x) = σ̃ Art_L(x) σ̃^{-1}` の**捩れ部分での言い換え**そのものである。

## ★何が新しく分かったか(★測定結果。楽観的に書かない)

| 段 | 内容 | 本ファイルの結果 |
|---|---|---|
| (a) | `Λ_n(S) ≅ μ_{p^n}(K̄)`(同型があること) | ★**証明した**(`nonempty_cyclotomeEquivRootsOfUnityClosure`) |
| (b) | `g ∈ Γ_F` の `μ_{p^n}(K̄)` への作用は `ζ ↦ ζ^{χ_n(g)}` | ★**証明した**(`rootsOfUnityConj_eq_pow`、★**定義そのもの**) |
| (c) | (a) の同型が `Γ_F`-同変に取れること | ★**残った**(`ArtinEquivariance`) |
| (d) | 壁の言い方が**同値**であること | ★**証明した**(下記 5 形が全部同値) |
| (e) | `g` は `L_S` の自己同型 `σ_g` を誘導する | ★**証明した**(`fixedFieldAut`) |
| (f) | `μ_{p^n}(L_S) = μ_{p^n}(K̄)`(`S ≤ muFixer` の下) | ★**証明した**(`rootsOfUnityFixedFieldEquiv`、同変) |
| (g) | 群論の共役 = `Gal(K̄/L_S)` の `σ_g`-半線型共役 | ★**証明した**(`fixedFieldGalCME_conjSubgroupCME`) |

## ★★壁の 5 つの同値な形(★どれを埋めても Proposition 1.1 が閉じる)

| 名前 | 主語 | 接ぎ先 |
|---|---|---|
| `CyclotomeConjIsCyclotomic`(Λ11) | `Λ_n(S)` の上で共役が `x^{χ_n(g)}` | —— |
| `ArtinEquivariance` | `Λ_n(S) ≅ μ_{p^n}(K̄)` が `Γ_F`-同変 | §4 |
| `ArtinEquivarianceFixedField` | `Λ_n(S) ≅ μ_{p^n}(L_S)` が `σ_g`-同変 | §5 |
| `ArtinEquivarianceGal` | `tors_{p^n}(Gal(K̄/L_S)^{ab}) ≅ μ_{p^n}(L_S)` が `τ ↦ gτg^{-1}` と同変 | §6 |
| ★**`ArtinEquivarianceLocalField`** | `tors_{p^n}(Γ_{L_S}^{ab}) ≅ μ_{p^n}(L_S)` が同変 | §7 ★**Λ9 に直結** |

★同値性は `artinEquivariance_iff_cyclotomeConj` /
`artinEquivarianceFixedField_iff_cyclotomeConj` /
`artinEquivarianceGal_iff_cyclotomeConj` /
`artinEquivarianceLocalField_iff_cyclotomeConj` の 4 本。

★(b) は `cyclotomicCharacter.spec`(mathlib)であって、こちら側には内容が無い。
★★**したがって壁は (c) だけであり、それは「Artin 写像が体の自己同型と同変」以外の
何ものでもない。** Λ11 が「指数そのものの決定」と書いた内容が、
本ファイルで「同型を同変に取れるか」という**1 つの選択の問題**に置き換わった。

## ★★なぜ `μ_{p^n}(L_S)` ではなく `μ_{p^n}(K̄)` で書くのか(★設計の要点)

Λ9 の到達点は `tors_{p^n}(Gal(L^{ab}/L)) ≅ μ_{p^n}(L)` であり、右辺は
**`L = L_S` の中の**冪根の群である。ところが `Γ_F` は `L` に作用しても
`μ_{p^n}(L)` の元を動かすので、右辺を `Γ_F`-加群として扱うには
`L` の中に閉じこもってはいけない。

★本ファイルは右辺を **`K̄ = F.closure` の中の `μ_{p^n}`** に取る。利点は 3 つ:

1. `Γ_F` の作用が**恒等的に**書ける(`g` を係数に当てるだけ、
   `MulEquiv.restrictRootsOfUnity`)。中間体を 1 枚も跨がないので **#59 に当たらない**。
2. その作用が `ζ ↦ ζ^{χ_n(g)}` であることが `cyclotomicCharacter.spec` から
   **そのまま**出る(`rootsOfUnityConj_eq_pow`)。
3. `S ≤ muFixer F (p^n)` のとき `μ_{p^n}(K̄) ⊆ L_S`
   (`rootsOfUnity_mem_fixedField`)なので、Λ9 の `μ_{p^n}(L_S)` との差は無い。
   ★**この 3 番目が「`μ_{p^n}(K̄)` で書いても弱くならない」ことの根拠**である。

## ★★★何が足りないか(★次の節点への申し送り)

`ArtinEquivariance` を埋めるには、`L := L_S` について

  「`Art_L : L^× → Gal(L^{ab}/L)` が `σ_g := g|_{L}` と同変」

が要る。★★**本ファイルは「同変性の左辺の翻訳」を全部済ませた** ——
残っているのは Artin 写像そのものの同変性 1 点である:

* `fixedFieldAut`(§5): `g ∈ Γ_F` は `S ⊴ Γ_F` のとき `L_S` の自己同型 `σ_g` を誘導する。
  ★`AlgEquiv.restrictNormalHom` を経由せず**直接**構成した(`Normal` インスタンスの
  探索も中間体の 2 層も要らない —— #59 の定型 (b))。
* `galConj`(§5): `Gal(K̄/L_S)` の上の `σ_g`-半線型共役 `τ ↦ g τ g^{-1}`。
* `fixedFieldGalCME_conjSubgroupCME`(§6): ★**群論の共役 `conjSubgroupCME g` は
  体の側では `galConj S g` に一致する** —— これが「群論側」と「体の側」を繋ぐ 1 本である。
* `cyclotomeEquiv_fixedFieldGalCME_conj`(§6)/`cyclotomeEquiv_absGalConjCME`(§7):
  その一致を円分子の上に持ち上げたもの。

★残るのは Lubin-Tate 側:`σ_g(π)` は別の素元であり、`σ_g` は形式群 `F_f` を
`F_{f^{σ_g}}` に、捩れ点 `μ_{f,m}` を `μ_{f^{σ_g},m}` に運ぶ。
★冪級数のレベルの道具は `Found/PGC/LubinTateEndoTwisted.lean`
(`LubinTateEndoTwisted_functional_equation`・`subst_LubinTateEndoTwisted_LubinTateEndoTwisted`)が
既に用意している。★**足りないのは「点で評価する層」**である
(同ファイルの申し送り)。素元の取り替えは
`LubinTateUniformizerIndependence.lean` の
`lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer` が**体の等式としては**
与えているが、★**写像の等式としては与えていない**(Cor 4.9 の具体化が要る)。

## ★退化の自己検査

1. **`Λ_n(S)` は位数ちょうど `p^n`**(Λ11 の `natCard_cyclotome_eq`)であり、
   `μ_{p^n}(K̄)` も位数ちょうど `p^n`(`natCard_rootsOfUnity_closure`)。
   ★どちらも自明群ではないので `ArtinEquivariance` は空虚でない。
2. **`S ≤ muFixer F (p^n)` は落とせない**。落とすと `Λ_n(S)` の位数が `p^n` を
   割るだけになり(Λ9)、`μ_{p^n}(K̄)` と同型でなくなる。
3. **`g ∈ S` の場合は既に証明済み**(Λ11 の
   `cyclotomeConj_eq_pow_cyclotomicCharacter_of_mem`)。★本ファイルが名指しする壁は
   `g ∉ S`、すなわち `Gal(L_S/F)` の上の話である。
4. **`K = ℚ_p`(`p` 奇)でも主張は成り立つ**。このとき `μ_p(ℚ_p) = 1` だが、
   `μ_p(K̄)` は位数 `p` であり `S ≤ muFixer` が `L_S ⊇ μ_p` を強制するので
   `L_S ≠ ℚ_p` になる。★`μ_{p^n}(K̄)` で書いたことでこの退化が消えている
   ——`μ_{p^n}(L_S)` で書くと `L_S` に依存して自明になりうる。
5. **`n = 0` でも成り立つ**(両辺とも自明群)。
6. **式の向きを変えていない**——`artinEquivariance_iff_cyclotomeConj` が
   Λ11 の壁との同値を与える。★言い換えで強くも弱くもなっていない。

## 逸脱の記録

* 本ファイルは新しい仮定を 1 つも置いていない。5 つの `Prop` はすべて
  `CyclotomeConjIsCyclotomic` と**同値**であり、Proposition 1.1 への配線は
  Λ11 のものをそのまま使う。★**通すために statement をねじ曲げていない。**
* 原典(pGC §1)の論拠は局所 Tate 双対性であって局所類体論ではない。この逸脱は
  Λ10/Λ11 と同じく `ResearchPaper/pgc-goal.md` に記録済み。

## ★測定の記録(★次の agent が再実行できる形で)

* ★**素元非依存性(`LubinTateUniformizerIndependence`)は消費しなかった。**
  本ファイルの射程は「同変性の左辺の翻訳」までで、Lubin-Tate に触れていない。
  素元の取り替えが要るのは、残った 1 点(Artin 写像そのものの同変性)の内側である。
* ★**`Found/PGC/LubinTateEndoTwisted.lean` も import していない。**
  同ファイルは冪級数の等式(`[θ]_{f,f′}` のねじれ版)を与えるが、本ファイルには
  冪級数が 1 つも現れないため接点が無かった。★消費するのは次の節点である。
* ★**#158(同名で書いて `already been declared` を出させる)の結果は 0 件。**
  `grep -rnE "^(theorem|lemma|noncomputable def|def) (fixedFieldAut|galConj|galConjAux|`
  `fixedFieldGalCME|rootsOfUnityConj|absGalConjCME|fixedField_map_mem|rootsOfUnityFixedField)"`
  `lean/ABC3/` も 0 件。★木に既存の在庫は無かった。
* ★mathlib からは `MulEquiv.restrictRootsOfUnity` /
  `rootsOfUnityEquivOfPrimitiveRoots` / `cyclotomicCharacter.spec` /
  `mulEquivOfCyclicCardEq` / `rootsOfUnity.isCyclic` を引いた
  (`grep -n "rootsOfUnity" .cache/mathlib-index.txt | grep -iE "equiv|map|aut|gal"` の 1 回)。
  ★**「mathlib に無い」と書いた箇所は本ファイルに 1 つも無い。**
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

universe u v

/-! ## 1. ★抽象核 —— 「同変な同型」と「冪で作用する」は同じこと

★この節に体・付値・分岐・位相の語彙は 1 つも出てこない。純粋な群論である。 -/

section AbstractCore

variable {Λ M : Type*} [Group Λ] [Group M]

/-- **★抽象核 1** —— 同変な同型があり、標的側で作用が `c` 乗なら、源でも `c` 乗。

`e (σ x) = τ (e x)` かつ `τ y = y ^ c` ならば `σ x = x ^ c`。
★`e` の単射性しか使わない。 -/
theorem eq_pow_of_equivariant (e : Λ ≃* M) (σ : Λ ≃* Λ) (τ : M ≃* M) (c : ℕ)
    (he : ∀ x, e (σ x) = τ (e x)) (hτ : ∀ y, τ y = y ^ c) (x : Λ) :
    σ x = x ^ c := by
  apply e.injective
  rw [he, hτ, ← map_pow]

/-- **★抽象核 1'(逆向き)** —— 両側で作用が同じ `c` 乗なら、**任意の**同型が同変。

★これが「同型の取り方の不定性は問題にならない」の中身である
(Λ10 の `pow_eq_pow_of_equivariant` と同じ現象の、こちら向きの言い方)。
★`e` については何も仮定していない —— 群同型でありさえすればよい。 -/
theorem equivariant_of_eq_pow (e : Λ ≃* M) (σ : Λ ≃* Λ) (τ : M ≃* M) (c : ℕ)
    (hσ : ∀ x, σ x = x ^ c) (hτ : ∀ y, τ y = y ^ c) (x : Λ) :
    e (σ x) = τ (e x) := by
  rw [hσ, hτ, map_pow]

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

/-! ## 2. `μ_m(K̄)` への `Γ_K` の作用

★中間体を 1 枚も跨がない —— `g` を係数に当てるだけ。 -/

/-- **`g : Γ_K` の `μ_m(K̄)` への作用**。

`g` は `K̄` の環同型なので単数群の同型を誘導し、`m` 乗根の群を保つ。
mathlib の `MulEquiv.restrictRootsOfUnity` に代入しただけ。

★これは**定義**であって定理ではない —— `Γ_K` は `K̄` に作用しているのだから、
`μ_m(K̄)` にも作用する。★**`L_S` の中の `μ_m` ではこの作用が書けない**
(`Γ_F` は `L_S` を動かす)ことが、本ファイルが `K̄` の中で書く理由である。 -/
noncomputable def rootsOfUnityConj (K : PAdicLocalField p) (m : ℕ) (g : K.absGal) :
    ↥(rootsOfUnity m K.closure) ≃* ↥(rootsOfUnity m K.closure) :=
  MulEquiv.restrictRootsOfUnity (g.toRingEquiv.toMulEquiv) m

@[simp] theorem rootsOfUnityConj_coe (K : PAdicLocalField p) (m : ℕ) (g : K.absGal)
    (x : ↥(rootsOfUnity m K.closure)) :
    (((rootsOfUnityConj K m g x : ↥(rootsOfUnity m K.closure)) : (K.closure)ˣ) : K.closure)
      = g ((x : (K.closure)ˣ) : K.closure) :=
  MulEquiv.restrictRootsOfUnity_coe_apply _ _

theorem rootsOfUnityConj_one (K : PAdicLocalField p) (m : ℕ)
    (x : ↥(rootsOfUnity m K.closure)) : rootsOfUnityConj K m 1 x = x := by
  refine Subtype.ext (Units.ext ?_)
  rw [rootsOfUnityConj_coe]
  rfl

theorem rootsOfUnityConj_mul (K : PAdicLocalField p) (m : ℕ) (g₁ g₂ : K.absGal)
    (x : ↥(rootsOfUnity m K.closure)) :
    rootsOfUnityConj K m (g₁ * g₂) x = rootsOfUnityConj K m g₁ (rootsOfUnityConj K m g₂ x) := by
  refine Subtype.ext (Units.ext ?_)
  rw [rootsOfUnityConj_coe, rootsOfUnityConj_coe, rootsOfUnityConj_coe]
  rfl

/-- `μ_m(K̄)` の元は `K̄` の中で `m` 乗して `1`。 -/
theorem rootsOfUnity_pow_eq_one (K : PAdicLocalField p) (m : ℕ)
    (x : ↥(rootsOfUnity m K.closure)) :
    (((x : (K.closure)ˣ) : K.closure)) ^ m = 1 := by
  have h : (x : (K.closure)ˣ) ^ m = 1 := (mem_rootsOfUnity m (x : (K.closure)ˣ)).mp x.2
  have h2 := congrArg (Units.val) h
  rwa [Units.val_pow_eq_pow_val, Units.val_one] at h2

/-- `μ_m(K̄)` の中の冪は `K̄` の中の冪。 -/
theorem rootsOfUnity_coe_pow (K : PAdicLocalField p) (m c : ℕ)
    (x : ↥(rootsOfUnity m K.closure)) :
    (((x ^ c : ↥(rootsOfUnity m K.closure)) : (K.closure)ˣ) : K.closure)
      = (((x : (K.closure)ˣ) : K.closure)) ^ c := by
  rw [SubmonoidClass.coe_pow, Units.val_pow_eq_pow_val]

/-- **★★`Γ_K` は `μ_{p^n}(K̄)` に `χ_{K,n}` 倍で作用する**。

★これは `cyclotomicCharacter.spec`(mathlib)そのもの —— **円分指標の定義**である。
こちら側に数学的内容は無い。

★★**壁の「右半分」がここで完全に消える。** 残るのは
「`Λ_n(S) ≅ μ_{p^n}(K̄)` を `Γ_F`-同変に取れるか」だけである。 -/
theorem rootsOfUnityConj_eq_pow (K : PAdicLocalField p) (n : ℕ) (g : K.absGal)
    (x : ↥(rootsOfUnity (p ^ n) K.closure)) :
    rootsOfUnityConj K (p ^ n) g x
      = x ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  refine Subtype.ext (Units.ext ?_)
  rw [rootsOfUnityConj_coe, rootsOfUnity_coe_pow]
  exact cyclotomicCharacter.spec (n := n) p g.toRingEquiv _ (rootsOfUnity_pow_eq_one K (p ^ n) x)

/-! ## 3. `μ_{p^n}(K̄)` は位数ちょうど `p^n` の巡回群 -/

theorem natCard_rootsOfUnity_closure (K : PAdicLocalField p) (n : ℕ) :
    Nat.card ↥(rootsOfUnity (p ^ n) K.closure) = p ^ n := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_closure' K (p ^ n)
  rw [Nat.card_eq_fintype_card, hζ.card_rootsOfUnity]

/-- **`Λ_n(S) ≅ μ_{p^n}(K̄)`**(同型があること)。

Λ11 の `nonempty_cyclotomeEquivZMod`(位数ちょうど `p^n` の巡回群)と、
`μ_{p^n}(K̄)` が同じく位数ちょうど `p^n` の巡回群であることの合成。
★**作用については何も言っていない** —— それが `ArtinEquivariance` である。 -/
theorem nonempty_cyclotomeEquivRootsOfUnityClosure (K : PAdicLocalField p) {n : ℕ}
    {S : Subgroup K.absGal} (hS : IsOpen (S : Set K.absGal))
    (hle : S ≤ muFixer K (p ^ n)) :
    Nonempty (↥(cyclotome ↥S (p ^ n)) ≃* ↥(rootsOfUnity (p ^ n) K.closure)) := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  haveI : IsCyclic ↥(rootsOfUnity (p ^ n) K.closure) := rootsOfUnity.isCyclic _ _
  obtain ⟨e⟩ := nonempty_cyclotomeEquivZMod K hS hle
  haveI : IsCyclic ↥(cyclotome ↥S (p ^ n)) :=
    (MulEquiv.isCyclic e).mpr (by infer_instance)
  refine ⟨mulEquivOfCyclicCardEq ?_⟩
  rw [natCard_cyclotome_eq K hS hle, natCard_rootsOfUnity_closure K n]

/-! ## 4. ★★★★★残った壁 —— Artin 写像の同変性 -/

def ArtinEquivariance.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★これだけが残った壁**(Λ11 の `CyclotomeConjIsCyclotomic` の同値な言い換え)。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

「`μ_{p^n} ⊆ L_S` なる開正規 `S ⊴ Γ_F` について、群論的円分子 `Λ_n(S) = tors_{p^n}(S^{ab})`
から `μ_{p^n}(K̄)` への **`Γ_F`-同変な**同型が存在する。」

★★これは Artin 写像の同変性 `Art_L(σ x) = σ̃ Art_L(x) σ̃^{-1}` の捩れ部分での言い換えである。
★同型の**存在**は既に証明済み(`nonempty_cyclotomeEquivRootsOfUnityClosure`)。
足りないのは**同変に取れること**だけである。

★退化の自己検査: 両辺とも位数ちょうど `p^n` の巡回群
(`natCard_cyclotome_eq` / `natCard_rootsOfUnity_closure`)なので空虚ではない。
`g ∈ S` の場合は両辺とも恒等になり、Λ11 の
`cyclotomeConj_eq_pow_cyclotomicCharacter_of_mem` が既に押さえている。 -/
def ArtinEquivariance (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal),
    IsOpen (S : Set F.absGal) → S ≤ muFixer F (p ^ n) →
    ∃ e : ↥(cyclotome ↥S (p ^ n)) ≃* ↥(rootsOfUnity (p ^ n) F.closure),
      ∀ (g : F.absGal) (x : ↥(cyclotome ↥S (p ^ n))),
        e (@cyclotomeConj F.absGal _ _ _ S hS (p ^ n) g x)
          = rootsOfUnityConj F (p ^ n) g (e x)

/-- **★★同変性 ⇒ Λ11 の壁**。抽象核 1 + `rootsOfUnityConj_eq_pow`。 -/
theorem cyclotomeConjIsCyclotomic_of_artinEquivariance
    (h : ArtinEquivariance p) : CyclotomeConjIsCyclotomic p := by
  intro F n S hS hopen hle g x
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  exact eq_pow_of_equivariant e _ (rootsOfUnityConj F (p ^ n) g) _ (he g)
    (rootsOfUnityConj_eq_pow F n g) x

/-- **★★Λ11 の壁 ⇒ 同変性**。抽象核 1' —— 作用が両側で `c` 乗なら**どの**同型も同変。

★これで「`μ_{p^n}(K̄)` で書き換えたことによって主張が強くなっていない」ことが確定する。 -/
theorem artinEquivariance_of_cyclotomeConjIsCyclotomic
    (h : CyclotomeConjIsCyclotomic p) : ArtinEquivariance p := by
  intro F n S hS hopen hle
  obtain ⟨e⟩ := nonempty_cyclotomeEquivRootsOfUnityClosure F hopen hle
  refine ⟨e, fun g x => ?_⟩
  exact equivariant_of_eq_pow e _ (rootsOfUnityConj F (p ^ n) g) _ (h F n S hS hopen hle g)
    (rootsOfUnityConj_eq_pow F n g) x

/-- **★★★2 つの言い方は同値**。 -/
theorem artinEquivariance_iff_cyclotomeConj :
    ArtinEquivariance p ↔ CyclotomeConjIsCyclotomic p :=
  ⟨cyclotomeConjIsCyclotomic_of_artinEquivariance,
    artinEquivariance_of_cyclotomeConjIsCyclotomic⟩

def cyclotomicCharacter_recoverable_of_artinEquivariance.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★★★★★★★★★[pGC] Proposition 1.1 の現在地**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

円分指標の群論的復元は、いまや **Artin 写像の同変性** 1 本だけに依存している。
☆★**2026-09-07 に埋まった**（`Found/PGC/ReciprocityDatumIndependence.lean::cyclotomicCharacter_recoverable_holds`、★仮定ゼロ）。★本ファイルに本物の `sorry` は 1 件も無い。 -/
theorem cyclotomicCharacter_recoverable_of_artinEquivariance
    (h : ArtinEquivariance p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicCharacter_recoverable_of_cyclotomeConj
    (cyclotomeConjIsCyclotomic_of_artinEquivariance h)

/-! ## 5. ★★「`σ` が体を運ぶ」段 —— `g ∈ Γ_F` から `σ_g : L_S ≃ L_S` へ

★この節の §5.1 に局所体・付値・分岐の語彙は 1 つも出てこない(一般の体拡大の話)。
★★**次の節点(Lubin-Tate 側)への入口はここである** —— Artin 写像の同変性を
述べるには、まず `g` が誘導する `L_S` の自己同型 `σ_g` を持たねばならない。 -/

section SigmaCarries

/-! ### 5.1 抽象核 —— 正規部分群の固定体は安定 -/

variable {k E : Type u} [Field k] [Field E] [Algebra k E]

/-- **★抽象核 2** —— `S ⊴ Gal(E/k)` なら `Gal(E/k)` の元は `fixedField S` を保つ。

★分岐・付値・局所体の語彙ゼロ。使うのは `S` の正規性だけである。
証明は `g⁻¹ f g ∈ S` の 1 行(正規性)と、それを `x` に当てて `g` で戻すこと。 -/
theorem fixedField_map_mem (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E) {x : E}
    (hx : x ∈ IntermediateField.fixedField S) :
    g x ∈ IntermediateField.fixedField S := by
  rw [IntermediateField.mem_fixedField_iff] at hx ⊢
  intro f hf
  have hmem : g⁻¹ * f * g ∈ S := by
    have := Subgroup.Normal.conj_mem ‹S.Normal› f hf g⁻¹
    simpa using this
  have h2 : g⁻¹ (f (g x)) = x := hx _ hmem
  have h3 := congrArg (fun z => g z) h2
  rw [show g (g⁻¹ (f (g x))) = (g * g⁻¹) (f (g x)) from rfl, mul_inv_cancel] at h3
  exact h3

/-- **★★抽象核 3 —— `σ_g`**。`S ⊴ Gal(E/k)` のとき `g` が誘導する `fixedField S` の自己同型。

★`AlgEquiv.restrictNormalHom` を**経由しない**。`Normal k ↥(fixedField S)` の
インスタンス探索も、中間体を 2 層跨ぐ `rfl` も要らない(#59 の定型 (b):
`Gal(E/k)` から 1 層で書く)。 -/
noncomputable def fixedFieldAut (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E) :
    ↥(IntermediateField.fixedField S) ≃ₐ[k] ↥(IntermediateField.fixedField S) where
  toFun x := ⟨g x, fixedField_map_mem S g x.2⟩
  invFun x := ⟨g⁻¹ x, fixedField_map_mem S g⁻¹ x.2⟩
  left_inv := by
    intro x
    refine Subtype.ext ?_
    show g⁻¹ (g (x : E)) = (x : E)
    rw [show g⁻¹ (g (x : E)) = (g⁻¹ * g) (x : E) from rfl, inv_mul_cancel]
    rfl
  right_inv := by
    intro x
    refine Subtype.ext ?_
    show g (g⁻¹ (x : E)) = (x : E)
    rw [show g (g⁻¹ (x : E)) = (g * g⁻¹) (x : E) from rfl, mul_inv_cancel]
    rfl
  map_mul' := fun x y => Subtype.ext (map_mul g _ _)
  map_add' := fun x y => Subtype.ext (map_add g _ _)
  commutes' := fun r => Subtype.ext (g.commutes r)

@[simp] theorem coe_fixedFieldAut (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E)
    (x : ↥(IntermediateField.fixedField S)) :
    ((fixedFieldAut S g x : ↥(IntermediateField.fixedField S)) : E) = g (x : E) := rfl

theorem fixedFieldAut_one (S : Subgroup (E ≃ₐ[k] E)) [S.Normal]
    (x : ↥(IntermediateField.fixedField S)) : fixedFieldAut S 1 x = x :=
  Subtype.ext rfl

theorem fixedFieldAut_mul (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g₁ g₂ : E ≃ₐ[k] E)
    (x : ↥(IntermediateField.fixedField S)) :
    fixedFieldAut S (g₁ * g₂) x = fixedFieldAut S g₁ (fixedFieldAut S g₂ x) :=
  Subtype.ext rfl

/-- ★**`g ∈ S` なら `σ_g` は恒等** —— 作用は `Gal(E/k) ⧸ S` を経由する。

★これが「壁の易しい半分」の体の側の姿である(群の側は Λ2 の
`cyclotomeConj_coe_self`)。 -/
theorem fixedFieldAut_of_mem (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] {g : E ≃ₐ[k] E}
    (hg : g ∈ S) (x : ↥(IntermediateField.fixedField S)) : fixedFieldAut S g x = x :=
  Subtype.ext ((IntermediateField.mem_fixedField_iff _ _).mp x.2 g hg)

/-! ### 5.1b ★`Gal(E/L_S)` の上の `σ_g`-半線型共役

★**Artin 写像の同変性が主語にしている作用そのもの**である:
`τ ↦ g τ g^{-1}` は `Gal(E/L_S)` の群自己同型で、`L_S` の上では `σ_g` として働く
(だから「`σ_g`-半線型」)。★この節も分岐・付値・局所体の語彙ゼロ。 -/

/-- `g τ g^{-1}` は `L_S` 上線型である —— `g^{-1}` が `L_S` を保つから。 -/
noncomputable def galConjAux (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E)
    (τ : E ≃ₐ[↥(IntermediateField.fixedField S)] E) :
    E ≃ₐ[↥(IntermediateField.fixedField S)] E where
  toFun x := g (τ (g⁻¹ x))
  invFun x := g (τ.symm (g⁻¹ x))
  left_inv := by
    intro x
    show g (τ.symm (g⁻¹ (g (τ (g⁻¹ x))))) = x
    rw [show g⁻¹ (g (τ (g⁻¹ x))) = (g⁻¹ * g) (τ (g⁻¹ x)) from rfl, inv_mul_cancel]
    show g (τ.symm (τ (g⁻¹ x))) = x
    rw [τ.symm_apply_apply, show g (g⁻¹ x) = (g * g⁻¹) x from rfl, mul_inv_cancel]
    rfl
  right_inv := by
    intro x
    show g (τ (g⁻¹ (g (τ.symm (g⁻¹ x))))) = x
    rw [show g⁻¹ (g (τ.symm (g⁻¹ x))) = (g⁻¹ * g) (τ.symm (g⁻¹ x)) from rfl, inv_mul_cancel]
    show g (τ (τ.symm (g⁻¹ x))) = x
    rw [τ.apply_symm_apply, show g (g⁻¹ x) = (g * g⁻¹) x from rfl, mul_inv_cancel]
    rfl
  map_mul' := by
    intro x y
    show g (τ (g⁻¹ (x * y))) = g (τ (g⁻¹ x)) * g (τ (g⁻¹ y))
    rw [map_mul, map_mul, map_mul]
  map_add' := by
    intro x y
    show g (τ (g⁻¹ (x + y))) = g (τ (g⁻¹ x)) + g (τ (g⁻¹ y))
    rw [map_add, map_add, map_add]
  commutes' := by
    intro r
    show g (τ (g⁻¹ (r : E))) = (r : E)
    have hmem : g⁻¹ (r : E) ∈ IntermediateField.fixedField S := fixedField_map_mem S g⁻¹ r.2
    have hfix : τ (g⁻¹ (r : E)) = g⁻¹ (r : E) := τ.commutes ⟨g⁻¹ (r : E), hmem⟩
    rw [hfix, show g (g⁻¹ (r : E)) = (g * g⁻¹) (r : E) from rfl, mul_inv_cancel]
    rfl

@[simp] theorem galConjAux_apply (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E)
    (τ : E ≃ₐ[↥(IntermediateField.fixedField S)] E) (x : E) :
    galConjAux S g τ x = g (τ (g⁻¹ x)) := rfl

/-- **★★抽象核 4 —— `σ_g`-半線型共役** `τ ↦ g τ g^{-1}` を `Gal(E/L_S)` の群自己同型として。

★★**Artin 写像の同変性 `Art(σ x) = g Art(x) g^{-1}` の右辺に現れる作用がこれ**である。 -/
noncomputable def galConj (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E) :
    (E ≃ₐ[↥(IntermediateField.fixedField S)] E)
      ≃* (E ≃ₐ[↥(IntermediateField.fixedField S)] E) where
  toFun := galConjAux S g
  invFun := galConjAux S g⁻¹
  left_inv := by
    intro τ
    refine AlgEquiv.ext (fun x => ?_)
    show g⁻¹ (g (τ (g⁻¹ ((g⁻¹)⁻¹ x)))) = τ x
    rw [inv_inv, show g⁻¹ (g (τ (g⁻¹ (g x)))) = (g⁻¹ * g) (τ (g⁻¹ (g x))) from rfl, inv_mul_cancel]
    show τ (g⁻¹ (g x)) = τ x
    rw [show g⁻¹ (g x) = (g⁻¹ * g) x from rfl, inv_mul_cancel]
    rfl
  right_inv := by
    intro τ
    refine AlgEquiv.ext (fun x => ?_)
    show g (g⁻¹ (τ ((g⁻¹)⁻¹ (g⁻¹ x)))) = τ x
    rw [inv_inv, show g (g⁻¹ (τ (g (g⁻¹ x)))) = (g * g⁻¹) (τ (g (g⁻¹ x))) from rfl, mul_inv_cancel]
    show τ (g (g⁻¹ x)) = τ x
    rw [show g (g⁻¹ x) = (g * g⁻¹) x from rfl, mul_inv_cancel]
    rfl
  map_mul' := by
    intro τ₁ τ₂
    refine AlgEquiv.ext (fun x => ?_)
    show g (τ₁ (τ₂ (g⁻¹ x))) = g (τ₁ (g⁻¹ (g (τ₂ (g⁻¹ x)))))
    rw [show g⁻¹ (g (τ₂ (g⁻¹ x))) = (g⁻¹ * g) (τ₂ (g⁻¹ x)) from rfl, inv_mul_cancel]
    rfl

@[simp] theorem galConj_apply (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E)
    (τ : E ≃ₐ[↥(IntermediateField.fixedField S)] E) (x : E) :
    galConj S g τ x = g (τ (g⁻¹ x)) := rfl

/-- ★**`galConj` は `L_S` の上で `σ_g` として働く** —— 「半線型」の中身。

★これが Artin 写像の同変性の主語である
「`Art_{L}(σ_g x) = g Art_L(x) g^{-1}`」の両辺を結ぶ関係式である。 -/
theorem galConj_comm_fixedFieldAut (S : Subgroup (E ≃ₐ[k] E)) [S.Normal] (g : E ≃ₐ[k] E)
    (τ : E ≃ₐ[↥(IntermediateField.fixedField S)] E)
    (x : ↥(IntermediateField.fixedField S)) :
    galConj S g τ ((fixedFieldAut S g x : ↥(IntermediateField.fixedField S)) : E)
      = g ((τ (x : E))) := by
  show g (τ (g⁻¹ (g (x : E)))) = g (τ (x : E))
  rw [show g⁻¹ (g (x : E)) = (g⁻¹ * g) (x : E) from rfl, inv_mul_cancel]
  rfl

end SigmaCarries

/-! ### 5.2 `μ_{p^n}(L_S) = μ_{p^n}(K̄)` —— `S ≤ muFixer` の下で -/

/-- `S ≤ muFixer K m` なら `L_S = fixedField S` の中に原始 `m` 乗根がある。

★Λ10 の `exists_isPrimitiveRoot_fixedField_of_le_muFixer` と同じものだが、
**開性を仮定しない**形に整えてある(`fixedFieldLocalField` を経由しないので
`IsOpen` が要らない)。 -/
theorem exists_isPrimitiveRoot_fixedField' (K : PAdicLocalField p) {m : ℕ} [NeZero m]
    {S : Subgroup K.absGal} (hle : S ≤ muFixer K m) :
    ∃ η : ↥(IntermediateField.fixedField S), IsPrimitiveRoot η m := by
  obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_closure' K m
  have hmem : ζ ∈ IntermediateField.fixedField S :=
    (le_muFixer_iff K m S).mp hle ζ hζ.pow_eq_one
  exact ⟨⟨ζ, hmem⟩, IsPrimitiveRoot.of_map_of_injective
    (f := (IntermediateField.fixedField S).val) hζ (fun _ _ h => Subtype.ext h)⟩

theorem primitiveRoots_fixedField_nonempty (K : PAdicLocalField p) {m : ℕ} [NeZero m]
    {S : Subgroup K.absGal} (hle : S ≤ muFixer K m) :
    (primitiveRoots m ↥(IntermediateField.fixedField S)).Nonempty := by
  obtain ⟨η, hη⟩ := exists_isPrimitiveRoot_fixedField' K hle
  exact ⟨η, (mem_primitiveRoots (Nat.pos_of_ne_zero (NeZero.ne m))).mpr hη⟩

/-- **★★`μ_m(L_S) ≅ μ_m(K̄)`**(`S ≤ muFixer K m` の下で)。

`L_S` は原始 `m` 乗根を含むので、包含 `L_S ↪ K̄` は `m` 乗根の群の**同型**を誘導する。
★**これが「`μ_m(K̄)` で書いても弱くならない」ことの根拠**(モジュール docstring の 3)。 -/
noncomputable def rootsOfUnityFixedFieldEquiv (K : PAdicLocalField p) {m : ℕ} [NeZero m]
    {S : Subgroup K.absGal} (hle : S ≤ muFixer K m) :
    ↥(rootsOfUnity m ↥(IntermediateField.fixedField S)) ≃* ↥(rootsOfUnity m K.closure) :=
  rootsOfUnityEquivOfPrimitiveRoots (f := (IntermediateField.fixedField S).val)
    (fun _ _ h => Subtype.ext h) (primitiveRoots_fixedField_nonempty K hle)

@[simp] theorem coe_rootsOfUnityFixedFieldEquiv (K : PAdicLocalField p) {m : ℕ} [NeZero m]
    {S : Subgroup K.absGal} (hle : S ≤ muFixer K m)
    (x : ↥(rootsOfUnity m ↥(IntermediateField.fixedField S))) :
    (((rootsOfUnityFixedFieldEquiv K hle x : ↥(rootsOfUnity m K.closure))
        : (K.closure)ˣ) : K.closure)
      = (((x : (↥(IntermediateField.fixedField S))ˣ) : ↥(IntermediateField.fixedField S))
          : K.closure) := by
  have h := rootsOfUnityEquivOfPrimitiveRoots_symm_apply
    (f := (IntermediateField.fixedField S).val) (fun _ _ h => Subtype.ext h)
    (primitiveRoots_fixedField_nonempty K hle) (rootsOfUnityFixedFieldEquiv K hle x)
  rw [show (rootsOfUnityEquivOfPrimitiveRoots (f := (IntermediateField.fixedField S).val)
      (fun _ _ h => Subtype.ext h) (primitiveRoots_fixedField_nonempty K hle)).symm
      (rootsOfUnityFixedFieldEquiv K hle x) = x from
    (rootsOfUnityFixedFieldEquiv K hle).symm_apply_apply x] at h
  exact h.symm

/-- **`σ_g` の `μ_m(L_S)` への作用**。 -/
noncomputable def rootsOfUnityFixedFieldConj (K : PAdicLocalField p) (m : ℕ)
    (S : Subgroup K.absGal) [S.Normal] (g : K.absGal) :
    ↥(rootsOfUnity m ↥(IntermediateField.fixedField S))
      ≃* ↥(rootsOfUnity m ↥(IntermediateField.fixedField S)) :=
  MulEquiv.restrictRootsOfUnity (fixedFieldAut S g).toRingEquiv.toMulEquiv m

@[simp] theorem rootsOfUnityFixedFieldConj_coe (K : PAdicLocalField p) (m : ℕ)
    (S : Subgroup K.absGal) [S.Normal] (g : K.absGal)
    (x : ↥(rootsOfUnity m ↥(IntermediateField.fixedField S))) :
    (((rootsOfUnityFixedFieldConj K m S g x
          : ↥(rootsOfUnity m ↥(IntermediateField.fixedField S)))
        : (↥(IntermediateField.fixedField S))ˣ) : ↥(IntermediateField.fixedField S))
      = fixedFieldAut S g (((x : (↥(IntermediateField.fixedField S))ˣ)
          : ↥(IntermediateField.fixedField S))) :=
  MulEquiv.restrictRootsOfUnity_coe_apply _ _

/-- **★★★`μ_m(L_S) ≅ μ_m(K̄)` は `Γ_K`-同変**。

`σ_g` は `g` の制限なのだから、包含が誘導する同型は当然 `g` の作用と可換である。
★**これで「`L_S` の中の `μ`」と「`K̄` の中の `μ`」の差が完全に消える。** -/
theorem rootsOfUnityFixedFieldEquiv_conj (K : PAdicLocalField p) {m : ℕ} [NeZero m]
    {S : Subgroup K.absGal} [S.Normal] (hle : S ≤ muFixer K m) (g : K.absGal)
    (x : ↥(rootsOfUnity m ↥(IntermediateField.fixedField S))) :
    rootsOfUnityFixedFieldEquiv K hle (rootsOfUnityFixedFieldConj K m S g x)
      = rootsOfUnityConj K m g (rootsOfUnityFixedFieldEquiv K hle x) := by
  refine Subtype.ext (Units.ext ?_)
  rw [coe_rootsOfUnityFixedFieldEquiv, rootsOfUnityConj_coe, coe_rootsOfUnityFixedFieldEquiv,
    rootsOfUnityFixedFieldConj_coe, coe_fixedFieldAut]

/-- ★**`σ_g` は `μ_{p^n}(L_S)` に `χ_{K,n}(g)` 倍で作用する**。

§2 の `rootsOfUnityConj_eq_pow`(`K̄` の中の話)を §5.2 の同型で `L_S` に降ろしたもの。
★**Λ9 の `hΨ` が要求している式の、`μ` 側の姿がこれである。** -/
theorem rootsOfUnityFixedFieldConj_eq_pow (K : PAdicLocalField p) (n : ℕ)
    {S : Subgroup K.absGal} [S.Normal] (hle : S ≤ muFixer K (p ^ n)) (g : K.absGal)
    (x : ↥(rootsOfUnity (p ^ n) ↥(IntermediateField.fixedField S))) :
    rootsOfUnityFixedFieldConj K (p ^ n) S g x
      = x ^ ((PadicInt.toZModPow n
          ((cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p]))).val) := by
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  refine (rootsOfUnityFixedFieldEquiv K hle).injective ?_
  rw [rootsOfUnityFixedFieldEquiv_conj K hle g x, rootsOfUnityConj_eq_pow, map_pow]

/-! ### 5.3 ★壁を `L_S` の中の `μ` で書き直す(★次の節点が使う形) -/

/-- **★★★壁の `L_S` 版** —— Λ9 の到達点
(`nonempty_cyclotomeSubgroupEquivRootsOfUnity : Λ_n(S) ≅ μ_{p^n}(L_S)`)に
**同変性を足しただけ**の形。

★★**次の節点はこの形を埋めればよい** —— Λ9 が構成した同型
`torsionEquivRootsOfUnity` が `σ_g` と同変であること、すなわち
Artin 写像 `Art_{L_S}` が `σ_g` と同変であることを示せばよい。 -/
def ArtinEquivarianceFixedField (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal),
    IsOpen (S : Set F.absGal) → S ≤ muFixer F (p ^ n) →
    ∃ e : ↥(cyclotome ↥S (p ^ n))
        ≃* ↥(rootsOfUnity (p ^ n) ↥(IntermediateField.fixedField S)),
      ∀ (g : F.absGal) (x : ↥(cyclotome ↥S (p ^ n))),
        e (@cyclotomeConj F.absGal _ _ _ S hS (p ^ n) g x)
          = @rootsOfUnityFixedFieldConj p _ F (p ^ n) S hS g (e x)

/-- **★`L_S` 版 ⇒ `K̄` 版**。§5.2 の同変な同型で写すだけ。 -/
theorem artinEquivariance_of_fixedField
    (h : ArtinEquivarianceFixedField p) : ArtinEquivariance p := by
  intro F n S hS hopen hle
  haveI := hS
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine ⟨e.trans (rootsOfUnityFixedFieldEquiv F hle), fun g x => ?_⟩
  show rootsOfUnityFixedFieldEquiv F hle (e _) = _
  rw [he g x]
  exact rootsOfUnityFixedFieldEquiv_conj F hle g (e x)

/-- **★`K̄` 版 ⇒ `L_S` 版**。逆向きも同じ同型で戻る。 -/
theorem artinEquivarianceFixedField_of_artinEquivariance
    (h : ArtinEquivariance p) : ArtinEquivarianceFixedField p := by
  intro F n S hS hopen hle
  haveI := hS
  haveI : NeZero (p ^ n) := ⟨pn_ne_zero p n⟩
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine ⟨e.trans (rootsOfUnityFixedFieldEquiv F hle).symm, fun g x => ?_⟩
  refine (rootsOfUnityFixedFieldEquiv F hle).injective ?_
  rw [rootsOfUnityFixedFieldEquiv_conj F hle g]
  show (rootsOfUnityFixedFieldEquiv F hle) ((rootsOfUnityFixedFieldEquiv F hle).symm (e _))
    = rootsOfUnityConj F (p ^ n) g
        ((rootsOfUnityFixedFieldEquiv F hle) ((rootsOfUnityFixedFieldEquiv F hle).symm (e x)))
  rw [MulEquiv.apply_symm_apply, MulEquiv.apply_symm_apply, he g x]

/-- **★★★3 つの言い方は同値**。★言い換えで主張を強くも弱くもしていない。 -/
theorem artinEquivarianceFixedField_iff_cyclotomeConj :
    ArtinEquivarianceFixedField p ↔ CyclotomeConjIsCyclotomic p :=
  ⟨fun h => cyclotomeConjIsCyclotomic_of_artinEquivariance (artinEquivariance_of_fixedField h),
    fun h => artinEquivarianceFixedField_of_artinEquivariance
      (artinEquivariance_of_cyclotomeConjIsCyclotomic h)⟩

/-! ## 6. ★★群論の共役 = 体の自己同型による共役

★★**本ファイルでいちばん実質のある節**。Λ2 が群論だけで定義した `cyclotomeConj`
(正規部分群への共役)が、体の側では `Gal(K̄/L_S)` の上の `σ_g`-半線型共役
`τ ↦ g τ g^{-1}` に**一致する**ことを示す。

★これで残った壁は、完全に古典的な形になる:

  「`Gal(K̄/L)^{ab}` の `p^n` 捩れから `μ_{p^n}(L)` への同型で、
   `τ ↦ g τ g^{-1}` を `ζ ↦ σ_g(ζ)` に移すものが存在する。」

★★**これが Artin 写像の同変性そのもの**である。 -/

/-- **`↥S ≃ₜ* Gal(K̄/L_S)`**(`S` 開)。

★`absGalFixedFieldCME`(`AdjoinFieldClosure.lean`)と違い、**代数閉包の同一視を
経由しない** —— `K̄` をそのまま `L_S` の代数閉包として使う。
★中間体の層が 1 枚しかないので #59 に当たらない(定型 (b))。 -/
noncomputable def fixedFieldGalCME (K : PAdicLocalField p) (S : Subgroup K.absGal)
    (hS : IsOpen (S : Set K.absGal)) :
    ContinuousMulEquiv ↥S (K.closure ≃ₐ[↥(IntermediateField.fixedField S)] K.closure) := by
  haveI := isGalois_closure K
  haveI := finiteDimensional_fixedField_of_isOpen K S hS
  have h2 : (IntermediateField.fixedField S).fixingSubgroup = S :=
    InfiniteGalois.fixingSubgroup_fixedField
      (⟨S, Subgroup.isClosed_of_isOpen S hS⟩ : ClosedSubgroup K.absGal)
  exact (subgroupCongrContinuousMulEquiv h2.symm).trans
    (fixingSubgroupContinuousMulEquiv (IntermediateField.fixedField S))

@[simp] theorem fixedFieldGalCME_apply (K : PAdicLocalField p) (S : Subgroup K.absGal)
    (hS : IsOpen (S : Set K.absGal)) (x : ↥S) (y : K.closure) :
    fixedFieldGalCME K S hS x y = (x : K.absGal) y := rfl

theorem fixedFieldGalCME_symm_apply (K : PAdicLocalField p) (S : Subgroup K.absGal)
    (hS : IsOpen (S : Set K.absGal))
    (τ : K.closure ≃ₐ[↥(IntermediateField.fixedField S)] K.closure) (y : K.closure) :
    ((((fixedFieldGalCME K S hS).symm τ : ↥S)) : K.absGal) y = τ y := by
  conv_rhs => rw [← ContinuousMulEquiv.apply_symm_apply (fixedFieldGalCME K S hS) τ]
  rfl

/-- **★★★群論の共役は `σ_g`-半線型共役である**。

`conjSubgroupCME g`(Λ2、`s ↦ g s g^{-1}` を `↥S` の位相群同型として)を
`fixedFieldGalCME` で `Gal(K̄/L_S)` に移すと、§5 の `galConj S g`(`τ ↦ g τ g^{-1}`)に
なる。★**両辺とも `K̄` の上の写像としては同じ合成なので、`AlgEquiv.ext` 1 行で済む。** -/
theorem fixedFieldGalCME_conjSubgroupCME (K : PAdicLocalField p) (S : Subgroup K.absGal)
    [S.Normal] (hS : IsOpen (S : Set K.absGal)) (g : K.absGal) (x : ↥S) :
    fixedFieldGalCME K S hS (conjSubgroupCME g x)
      = galConj S g (fixedFieldGalCME K S hS x) := by
  refine AlgEquiv.ext (fun y => ?_)
  rw [fixedFieldGalCME_apply, conjSubgroupCME_coe, galConj_apply, fixedFieldGalCME_apply]
  rfl

/-- **`Gal(K̄/L_S)` の上の共役を位相群同型として**(連続性は `conjSubgroupCME` から移送)。 -/
noncomputable def galConjCME (K : PAdicLocalField p) (S : Subgroup K.absGal) [S.Normal]
    (hS : IsOpen (S : Set K.absGal)) (g : K.absGal) :
    ContinuousMulEquiv (K.closure ≃ₐ[↥(IntermediateField.fixedField S)] K.closure)
      (K.closure ≃ₐ[↥(IntermediateField.fixedField S)] K.closure) :=
  ((fixedFieldGalCME K S hS).symm.trans (conjSubgroupCME (H := S) g)).trans
    (fixedFieldGalCME K S hS)

/-- ★`galConjCME` は(連続性を忘れれば)`galConj` そのもの。 -/
theorem galConjCME_eq (K : PAdicLocalField p) (S : Subgroup K.absGal) [S.Normal]
    (hS : IsOpen (S : Set K.absGal)) (g : K.absGal)
    (τ : K.closure ≃ₐ[↥(IntermediateField.fixedField S)] K.closure) :
    galConjCME K S hS g τ = galConj S g τ := by
  show fixedFieldGalCME K S hS
      (conjSubgroupCME (H := S) g ((fixedFieldGalCME K S hS).symm τ)) = _
  rw [fixedFieldGalCME_conjSubgroupCME, ContinuousMulEquiv.apply_symm_apply]

/-- **★★★★★円分子の上でも同じ** —— `Λ_n(S)` への `g` の共役作用は、
`Gal(K̄/L_S)` の円分子の上の `galConj` 作用に一致する。

★★**これが本ファイルの「σ が運ぶ」段の到達点**である:
壁は「群論的な `Λ_n(S)`」の話から「`Gal(K̄/L_S)^{ab}` の捩れ」の話に完全に翻訳された。 -/
theorem cyclotomeEquiv_fixedFieldGalCME_conj (K : PAdicLocalField p) (m : ℕ)
    (S : Subgroup K.absGal) [S.Normal] (hS : IsOpen (S : Set K.absGal)) (g : K.absGal)
    (x : ↥(cyclotome ↥S m)) :
    cyclotomeEquiv (fixedFieldGalCME K S hS) m (cyclotomeConj S m g x)
      = cyclotomeEquiv (galConjCME K S hS g) m
          (cyclotomeEquiv (fixedFieldGalCME K S hS) m x) := by
  rw [cyclotomeConj, cyclotomeEquiv_trans, cyclotomeEquiv_trans]
  refine cyclotomeEquiv_congr (fun a => ?_) m x
  show fixedFieldGalCME K S hS (conjSubgroupCME (H := S) g a)
    = fixedFieldGalCME K S hS
        (conjSubgroupCME (H := S) g ((fixedFieldGalCME K S hS).symm
          (fixedFieldGalCME K S hS a)))
  rw [ContinuousMulEquiv.symm_apply_apply]

def ArtinEquivarianceGal.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★残った壁の、完全に古典的な形**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

「`μ_{p^n} ⊆ L` なる有限次 Galois 拡大 `L/F`(`L = L_S`)について、
`tors_{p^n}(Gal(K̄/L)^{ab})` から `μ_{p^n}(L)` への同型で、
`τ ↦ g τ g^{-1}` を `ζ ↦ σ_g(ζ)` に移すものが存在する。」

★★これが **Artin 写像の同変性** `Art_L(σ_g ζ) = g · Art_L(ζ) · g^{-1}` の
捩れ部分での言い換えである。★群論の語彙(`cyclotomeConj`・`TopologicalAbelianization` の
共役)は左辺から消え、体と Galois 群だけが残っている。 -/
def ArtinEquivarianceGal (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal)
    (hopen : IsOpen (S : Set F.absGal)), S ≤ muFixer F (p ^ n) →
    ∃ e : ↥(cyclotome (F.closure ≃ₐ[↥(IntermediateField.fixedField S)] F.closure) (p ^ n))
        ≃* ↥(rootsOfUnity (p ^ n) ↥(IntermediateField.fixedField S)),
      ∀ (g : F.absGal)
        (x : ↥(cyclotome (F.closure ≃ₐ[↥(IntermediateField.fixedField S)] F.closure) (p ^ n))),
        e (cyclotomeEquiv (@galConjCME p _ F S hS hopen g) (p ^ n) x)
          = @rootsOfUnityFixedFieldConj p _ F (p ^ n) S hS g (e x)

/-- **★Galois 版 ⇒ `L_S` 版**。§6 の翻訳を挟むだけ。 -/
theorem artinEquivarianceFixedField_of_gal
    (h : ArtinEquivarianceGal p) : ArtinEquivarianceFixedField p := by
  intro F n S hS hopen hle
  haveI := hS
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine ⟨(cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).trans e, fun g x => ?_⟩
  show e (cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n) _) = _
  rw [cyclotomeEquiv_fixedFieldGalCME_conj, he g]
  rfl

/-- **★`L_S` 版 ⇒ Galois 版**。逆向きも同じ翻訳で戻る。 -/
theorem artinEquivarianceGal_of_fixedField
    (h : ArtinEquivarianceFixedField p) : ArtinEquivarianceGal p := by
  intro F n S hS hopen hle
  haveI := hS
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine ⟨(cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).symm.trans e, fun g x => ?_⟩
  have hx : (cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).symm
      (cyclotomeEquiv (galConjCME F S hopen g) (p ^ n) x)
      = cyclotomeConj S (p ^ n) g
          ((cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).symm x) := by
    refine (cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).injective ?_
    rw [MulEquiv.apply_symm_apply,
      cyclotomeEquiv_fixedFieldGalCME_conj F (p ^ n) S hopen g
        ((cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).symm x),
      MulEquiv.apply_symm_apply]
  show e ((cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).symm
      (cyclotomeEquiv (galConjCME F S hopen g) (p ^ n) x))
    = rootsOfUnityFixedFieldConj F (p ^ n) S g
        (e ((cyclotomeEquiv (fixedFieldGalCME F S hopen) (p ^ n)).symm x))
  rw [hx, he g]

/-- **★★★4 つの言い方は同値**。★言い換えで主張を強くも弱くもしていない。 -/
theorem artinEquivarianceGal_iff_cyclotomeConj :
    ArtinEquivarianceGal p ↔ CyclotomeConjIsCyclotomic p :=
  ⟨fun h => artinEquivarianceFixedField_iff_cyclotomeConj.mp
      (artinEquivarianceFixedField_of_gal h),
    fun h => artinEquivarianceGal_of_fixedField
      (artinEquivarianceFixedField_iff_cyclotomeConj.mpr h)⟩

/-- **★★★★★★★★★★★★★★★★[pGC] Proposition 1.1 —— Galois 版の同変性から**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

☆★**2026-09-07 に埋まった**（`Found/PGC/ReciprocityDatumIndependence.lean::cyclotomicCharacter_recoverable_holds`、★仮定ゼロ）。★本ファイルに本物の `sorry` は 1 件も無い。 -/
theorem cyclotomicCharacter_recoverable_of_artinEquivarianceGal
    (h : ArtinEquivarianceGal p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicCharacter_recoverable_of_cyclotomeConj (artinEquivarianceGal_iff_cyclotomeConj.mp h)

/-! ## 7. ★★Λ9 が直接消費できる形(`fixedFieldLocalField` 版)

★★**次の節点が実際に手を付けるのはこの形である。** Λ9 の到達点
(`exists_abelianGalTorsion_equiv_rootsOfUnity`)と Λ11 の
`topAbelianizationEquivAbelianGal` は `PAdicLocalField` を引数に取るので、
`L := fixedFieldLocalField F S hopen` の言葉で書いておく方が接ぎやすい。 -/

/-- **`Γ_{L_S}` の上の共役作用**(`absGalFixedFieldCME` で移送したもの)。 -/
noncomputable def absGalConjCME (K : PAdicLocalField p) (S : Subgroup K.absGal) [S.Normal]
    (hS : IsOpen (S : Set K.absGal)) (g : K.absGal) :
    ContinuousMulEquiv (fixedFieldLocalField K S hS).absGal
      (fixedFieldLocalField K S hS).absGal :=
  ((absGalFixedFieldCME K S hS).trans (conjSubgroupCME (H := S) g)).trans
    (absGalFixedFieldCME K S hS).symm

/-- `absGalFixedFieldCME` は `absGalConjCME` を `cyclotomeConj` に移す。 -/
theorem cyclotomeEquiv_absGalConjCME (K : PAdicLocalField p) (m : ℕ) (S : Subgroup K.absGal)
    [S.Normal] (hS : IsOpen (S : Set K.absGal)) (g : K.absGal)
    (x : ↥(cyclotome (fixedFieldLocalField K S hS).absGal m)) :
    cyclotomeEquiv (absGalFixedFieldCME K S hS) m
        (cyclotomeEquiv (absGalConjCME K S hS g) m x)
      = cyclotomeConj S m g (cyclotomeEquiv (absGalFixedFieldCME K S hS) m x) := by
  rw [cyclotomeEquiv_trans, cyclotomeConj, cyclotomeEquiv_trans]
  refine cyclotomeEquiv_congr (fun a => ?_) m x
  show absGalFixedFieldCME K S hS ((absGalFixedFieldCME K S hS).symm
      (conjSubgroupCME (H := S) g (absGalFixedFieldCME K S hS a)))
    = conjSubgroupCME (H := S) g (absGalFixedFieldCME K S hS a)
  rw [ContinuousMulEquiv.apply_symm_apply]

def ArtinEquivarianceLocalField.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★残った壁の、Λ9 に接ぎやすい形**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

「`L := L_S` について、`tors_{p^n}(Γ_L^{ab})` から `μ_{p^n}(L)` への同型で、
`Γ_F` の共役作用を `σ_g` の作用に移すものが存在する。」

★Λ11 の `nonempty_cyclotomeSubgroupEquivRootsOfUnity`
(`Λ_n(S) ≅ μ_{p^n}(L_S)`、**作用なし**)に、**同変性だけを足した形**である。
★★**次の節点はこれを埋めればよい** —— Λ9 が構成した同型
`torsionEquivRootsOfUnity`(Artin 写像の捩れ部分)が `σ_g` と同変であることを示せばよい。 -/
def ArtinEquivarianceLocalField (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) (hS : S.Normal)
    (hopen : IsOpen (S : Set F.absGal)), S ≤ muFixer F (p ^ n) →
    ∃ e : ↥(cyclotome (fixedFieldLocalField F S hopen).absGal (p ^ n))
        ≃* ↥(rootsOfUnity (p ^ n) (fixedFieldLocalField F S hopen).carrier),
      ∀ (g : F.absGal)
        (x : ↥(cyclotome (fixedFieldLocalField F S hopen).absGal (p ^ n))),
        e (cyclotomeEquiv (@absGalConjCME p _ F S hS hopen g) (p ^ n) x)
          = @rootsOfUnityFixedFieldConj p _ F (p ^ n) S hS g (e x)

/-- **★`fixedFieldLocalField` 版 ⇒ `L_S` 版**。 -/
theorem artinEquivarianceFixedField_of_localField
    (h : ArtinEquivarianceLocalField p) : ArtinEquivarianceFixedField p := by
  intro F n S hS hopen hle
  haveI := hS
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine ⟨(cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).symm.trans e, fun g x => ?_⟩
  have hx : (cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).symm
      (cyclotomeConj S (p ^ n) g x)
      = cyclotomeEquiv (absGalConjCME F S hopen g) (p ^ n)
          ((cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).symm x) := by
    refine (cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).injective ?_
    rw [MulEquiv.apply_symm_apply,
      cyclotomeEquiv_absGalConjCME F (p ^ n) S hopen g
        ((cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).symm x),
      MulEquiv.apply_symm_apply]
  show e ((cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).symm
      (cyclotomeConj S (p ^ n) g x))
    = rootsOfUnityFixedFieldConj F (p ^ n) S g
        (e ((cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).symm x))
  rw [hx, he g]

/-- **★`L_S` 版 ⇒ `fixedFieldLocalField` 版**。 -/
theorem artinEquivarianceLocalField_of_fixedField
    (h : ArtinEquivarianceFixedField p) : ArtinEquivarianceLocalField p := by
  intro F n S hS hopen hle
  haveI := hS
  obtain ⟨e, he⟩ := h F n S hS hopen hle
  refine ⟨(cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)).trans e, fun g x => ?_⟩
  show e (cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n)
      (cyclotomeEquiv (absGalConjCME F S hopen g) (p ^ n) x))
    = rootsOfUnityFixedFieldConj F (p ^ n) S g
        (e (cyclotomeEquiv (absGalFixedFieldCME F S hopen) (p ^ n) x))
  rw [cyclotomeEquiv_absGalConjCME, he g]

/-- **★★★5 つの言い方はすべて同値**。★言い換えで主張を強くも弱くもしていない。 -/
theorem artinEquivarianceLocalField_iff_cyclotomeConj :
    ArtinEquivarianceLocalField p ↔ CyclotomeConjIsCyclotomic p :=
  ⟨fun h => artinEquivarianceFixedField_iff_cyclotomeConj.mp
      (artinEquivarianceFixedField_of_localField h),
    fun h => artinEquivarianceLocalField_of_fixedField
      (artinEquivarianceFixedField_iff_cyclotomeConj.mpr h)⟩

def cyclotomicCharacter_recoverable_of_artinEquivarianceLocalField.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★★★★★★★★★★★★★★★★[pGC] Proposition 1.1 の到達点**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

★★**Proposition 1.1 は「Λ9 の同型 `tors_{p^n}(Γ_L^{ab}) ≅ μ_{p^n}(L)` が
`σ_g` と同変に取れる」1 点だけに依存している。**
☆★**2026-09-07 に埋まった**（`Found/PGC/ReciprocityDatumIndependence.lean::cyclotomicCharacter_recoverable_holds`、★仮定ゼロ）。★本ファイルに本物の `sorry` は 1 件も無い。 -/
theorem cyclotomicCharacter_recoverable_of_artinEquivarianceLocalField
    (h : ArtinEquivarianceLocalField p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicCharacter_recoverable_of_cyclotomeConj
    (artinEquivarianceLocalField_iff_cyclotomeConj.mp h)

/-! ## 8. ★退化の自己検査 -/

/-- ★**`g ∈ S` の場合は既に成り立っている**(壁の易しい半分、`μ` 側の姿)。

`g ∈ S` なら `σ_g` は恒等(`fixedFieldAut_of_mem`)であり、
`χ_n(g) = 1`(Λ11 の `toZModPow_cyclotomicCharacter_eq_one_of_mem_muFixer`)。
★両辺が**独立の理由で**恒等になる。 -/
theorem rootsOfUnityFixedFieldConj_of_mem (K : PAdicLocalField p) (m : ℕ)
    {S : Subgroup K.absGal} [S.Normal] {g : K.absGal} (hg : g ∈ S)
    (x : ↥(rootsOfUnity m ↥(IntermediateField.fixedField S))) :
    rootsOfUnityFixedFieldConj K m S g x = x := by
  refine Subtype.ext (Units.ext ?_)
  rw [rootsOfUnityFixedFieldConj_coe, fixedFieldAut_of_mem S hg]

/-- ★**円分子と `μ_{p^n}(K̄)` はどちらも位数ちょうど `p^n`** ——
`ArtinEquivariance` が空虚に真な形ではないことの確認。 -/
theorem natCard_eq_of_le_muFixer (K : PAdicLocalField p) {n : ℕ}
    {S : Subgroup K.absGal} (hS : IsOpen (S : Set K.absGal))
    (hle : S ≤ muFixer K (p ^ n)) :
    Nat.card ↥(cyclotome ↥S (p ^ n)) = Nat.card ↥(rootsOfUnity (p ^ n) K.closure) := by
  rw [natCard_cyclotome_eq K hS hle, natCard_rootsOfUnity_closure K n]

end ABC3.Found.PGC
