import ABC3.Found.PGC.AbsGalRamificationFiltration
import ABC3.Found.PGC.InertiaIdentification
import ABC3.Found.PGC.UnramifiedClosureRoots

/-!
# 不分岐部分の除去と惰性群の全射性(Y19b + Y19c)

`Found/PGC/AbsGalRamificationFiltration.lean`(Y19)は
`RamificationFiltration p` を「各有限段の段データ `StageFiltration K.absGal`」に還元し、
その入口で **2 つの数学が要る**と名指しした:

* (a) 下付き分岐群 `G_i`(`i ≥ 0`)が**不分岐底変換 `L/K ↝ L/L₀` で変わらない**こと、
* (b) `L ⊆ L′` のとき `Gal(L′/L′₀) ↠ Gal(L/L₀)` が全射で、その核が
  Yoshida Corollary 6.13(i) の `H` に一致すること。

本ファイルはこの 2 つを扱う。`L₀ := L ⊓ K^ur`(`K^ur` は `K` の**最大不分岐拡大**
であって代数閉包ではない)。

## ★★本ファイルが原典に対応する項目を持たない理由

Yoshida は §6.1 で設定そのものを「完全分岐」と置いている:

原文 (Yoshida08 p.14):
> 6.1. Ramification groups. Let K′/K be a finite totally ramified Galois extension of local
> fields, and set G := Gal(K′/K). For a uniformizer π of K′, we have O[scr]_K′ = O[scr][π]
> by Lemma 5.11. We write v := v_K′ and q = |O[scr]/p[frak]| = |O[scr]_K′/p[frak]_K′|.

すなわち「一般の有限 Galois 拡大から完全分岐の設定へ落とす」段を原典は**書いていない**
(設定に畳んでいる)。本ファイルはその段を立てるための基盤ノードであり、
**原典に対応する主張が無いので `.src` を持たない**。

## 何が入ったか / 何が入っていないか(★正直な線引き)

**入ったもの**

* §1 **抽象核**(純イデアル論 + 純群論。分岐・付値・Galois の語彙が 1 つも出てこない):
  `B` を固定したまま群を `G ⊇ H` と縮めても、`𝔪.inertia G ≤ H` でありさえすれば
  `G_n(H)` は `H.subtype` で `G_n(G)` の上に**同型に写る**
  (`map_subtype_lowerRamificationGroup_of_inertia_le` / `lowerRamificationGroupMulEquiv`)。
  ★これが「不分岐底変換で `G_i` が変わらない」の中身である——底が `K` でも `L₀` でも
  整数環 `B = 𝒪_L` は**同じ環**で、変わるのは群だけだから。
* §2 **抽象核**(純群論): 全射 `f : Γ →* G` と固定部分群 `I ≤ Γ` について、
  `N′ ≤ ker f` なら `f(I ⊔ N′) = f(I)`(`map_sup_of_le_ker`)、
  その制限は全射で核は `ker f`(`quotientEquivMap`)、そして
  `(I ⊔ N′) · N = I ⊔ N`(`coe_sup_mul_coe_eq`)——これが Y19 の
  `StageFiltration.compat` の形そのものである。
* §3 **具体層**(§1 の橋): `G_0 = ker(Gal(K(x)/K) → Gal(k_{K(x)}/k))`
  (`lowerRamificationGroupAdjoin_zero_eq_ker_residueGalHom`)。
  ★`Ideal.inertia` による `G_0` の定義が**古典的な惰性群(剰余体への還元の核)**に
  一致することの確認。
* §4 **具体層**(§2): `inertiaGal K L := (絶対惰性群 I_K の Gal(L/K) への像)` を置き、
  - `L₀ = L ⊓ K^ur` が `inertiaGal K L` の**固定体**であること
    (`lift_fixedField_inertiaGal`。したがって `inertiaGal K L = Gal(L/L₀)`)、
  - `Γ_K` への引き戻しが `I_K ⊔ L.fixingSubgroup` であること(`comap_inertiaGal`)、
  - `L ⊆ L′` に対する全射性(`map_restrictNormalHom_sup_fixingSubgroup`)と
    その核(`quotientEquivMap` の具体化)。
* §5 **`Gal(L/L₀) ≤ G_0(L/K)` と、その帰結「`L/L₀` は完全分岐」**:
  - `inertiaGal_le_lowerRamificationGroupAdjoin_zero`(★Teichmüller 表現を使う。
    本ファイルで一番中身のある定理)、
  - `lowerRamificationGroup_inertiaGal_zero_eq_top`
    (`G_0(L/L₀) = ⊤`。★★**これが「不分岐部分の除去」の到達点** ——
    一般の有限 Galois `L/K` を、底を `L₀ = L ⊓ K^ur` へ移すことで
    **Yoshida §6.1 の設定(完全分岐)に落とせる**)、
  - 整合性検査: `L/K` が完全分岐なら `inertiaGal K L = ⊤ = G_0`
    (`inertiaGal_eq_lowerRamificationGroupAdjoin_zero_of_isTotallyRamified`)。
* §6 `I_K` から作った段データ `inertiaStageFiltration` とその極限
  (★**原典の `Γ_K^v` ではない**。下の警告を読むこと)。

**入っていないもの(★新ノード)**

★★**逆向きの包含 `G_0(L/K) ≤ inertiaGal K L`** は**埋まっていない**
(したがって `G_0 = Gal(L/L₀)` の等式も未証明)。★方向を取り違えないこと:
本ファイルにあるのは `Gal(L/L₀) ≤ G_0` である。

逆向きに必要なのは「還元射 `Gal(L/K) → Gal(k_L/k)` の像が `[L₀ : K]` 以上」、
すなわち古典的には**還元射の全射性**である。在庫の `residueGalHom_bijective` は
`IsUnramifiedAdjoin` を仮定しており一般の `L` には当たらない。詰めるには
`residueFieldHom`(`TotallyRamified.lean:644`、`K(z) ≤ K(x)` に対する剰余体の射)を
`Gal` の作用と両立させる補題が要る——中間体 2 層をまたぐので独立のノードにする。

★そのため §1 の具体層は仮定 `G_0 ≤ H` を**仮定のまま**受け取る形にしてある
(`lowerRamificationGroupAdjoin_map_subtype`)。

## ★本体の見当は当たったか

本体の見当は「`G_i ⊆ G_0 = 惰性群 = Gal(L/L₀)` なので一致する」だった。
実測の結果:

* 「`G_i ⊆ G_0` だから、`G_0 ≤ H` さえあれば `G_i(H) = G_i(G)`」——**当たり**。
  しかも `Ideal.inertia` で定義してあるおかげで `subgroupOf` が `rfl` になり、
  抽象核は 3 行で済んだ。
* 「`G_0 = 惰性群`」——**当たり**(§3。`G_0 = ker(residueGalHom)`)。
* 「`惰性群 = Gal(L/L₀)`」——★**半分だけ当たり**。`Gal(L/L₀) ≤ G_0` は §5 で証明した
  (Teichmüller)。逆向きは未証明である。
* ★**ただし、消費側が本当に要るのは半分のほうだった**——`L/L₀` が完全分岐であること
  (`G_0(L/L₀) = ⊤`)は `Gal(L/L₀) ≤ G_0` だけから出る。等式は要らない。
  ★本体の段取り「(a) `G_i` が不分岐底変換で不変」は**必要より強い要求**だった。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. §1 の主張を「`L₀ := L ⊓ K^ur`」ではなく「`𝔪.inertia G ≤ H` なる**任意の**部分群 `H`」で
   立てた。原典(が書いていない段)より**仮定が弱い**。`H = Gal(L/L₀)` を代入すれば
   目的の形になるが、その代入に必要な `G_0 ≤ Gal(L/L₀)` は本ファイルの外である。
2. §2 の全射性を「有限段の群 `Gal(L′/L′₀) → Gal(L/L₀)`」ではなく
   **`Γ_K` の部分群の言葉**(`I_K ⊔ N′ ↠ inertiaGal K L`)で書いた。
   Y19 の `StageFiltration` が `Γ` の部分群だけで書かれているので、
   中間体の塔(`↥L` 上の `↥L′`)を作らずに済む——`lean-idioms.md` #59 を回避するための
   設計判断であり、消費側にとっては同値である。
3. §6 の `inertiaStageFiltration` は `v ≥ 0` で `Γ^v = I_K` を与える。
   ★**原典の `Γ_K^v` ではない**(真の `Γ_K^v` は `v > 0` で野性惰性群に含まれるので
   もっと小さい)。Y19 の `trivialStageFiltration`(`v > 0` で `⊥`)が
   下からの退化 witness だったのに対し、これは**上からの近似**である。
   ★G2 の非空虚 witness にも `Skeleton/PGC/Section2.lean` の入力にも使ってはならない。

## 退化の自己検査

* ★★**`i ≥ 0` を落とすと §1 は偽**。`lowerRamificationGroup B G n = (𝔪^{n+1}).inertia G` は
  `n : ℕ` なので指数は常に `≥ 1` だが、「`G_{-1}`」に当たる `(𝔪^0).inertia G = ⊤.inertia G`
  は `inertia_top_eq_top` により `⊤` であり、`H` に含まれない。
  ⇒ 指数を `n` ではなく `n + 1` にしてある(= `i ≥ 0`)ことが本質的。
* ★**`i = 0` で成り立つか**: 成り立つ。`G_0 = 𝔪.inertia G`
  (`lowerRamificationGroup_zero_eq_inertia`)なので仮定 `𝔪.inertia G ≤ H` がそのまま
  `G_0 ≤ H` であり、§1 の結論は `n = 0` でも正しい。
* ★**`L/K` が Galois でないと `L₀` が定まらない**: §4 は `[Normal K.carrier ↥L]` を
  要求している(`AlgEquiv.restrictNormalHom` がそれ無しでは書けない)。
* ★**`K^ur` は最大不分岐拡大**であって代数閉包ではない。`unramifiedClosure K` は
  `IntermediateField K.carrier K.closure` であり `K.closure` とは別物。
  実際 `inertiaGal K L = ⊥` になるのは `L ≤ K^ur` のときだけである。
* ★**§2 の核は本当に Corollary 6.13(i) の `H` か**: `Γ_K` の言葉では
  `ker(Γ_K ↠ Gal(L/K)) = L.fixingSubgroup` なので、`I_K ⊔ L′.fixingSubgroup` から
  `inertiaGal K L` への全射の核は `(I_K ⊔ L′.fixingSubgroup) ⊓ L.fixingSubgroup`
  である(`quotientEquivMap`)。有限段に降ろすとこれは
  `Gal(L′/L′₀) ⊓ Gal(L′/L) = Gal(L′/ L·L′₀)` であり、
  `Gal(L′/L′₀) / Gal(L′/L·L′₀) ≅ Gal(L·L′₀/L′₀) ≅ Gal(L / L ⊓ L′₀) = Gal(L/L₀)`
  ——Corollary 6.13(i) の `G ⊵ H` の `H` に一致する。
  ★ただし最後の 2 つの同型(合成体の Galois 群)は**本ファイルでは証明していない**。
  `Γ_K` の言葉での形(`map_restrictNormalHom_sup_fixingSubgroup`)だけを出してある。
-/

namespace ABC3.Found.PGC

open scoped Pointwise

/-! ## §1 抽象核 —— 群を縮めても下付き分岐群は変わらない

★★分岐・付値・Galois の語彙が 1 つも出てこない。可換環 `B`(局所環)と
それに環として作用する群 `G` だけで書けている。 -/

/-- ★退化検査用: `⊤` の inertia は `⊤`。

原典の添字で言えば「`G_{-1}`」に当たるもので、これは惰性群より**大きい**。
⇒ §1 の主張は `i ≥ 0`(= 指数 `n + 1 ≥ 1`)を落とすと偽になる。 -/
theorem inertia_top_eq_top (B : Type*) [CommRing B] (G : Type*) [Group G]
    [MulSemiringAction G B] : (⊤ : Ideal B).inertia G = ⊤ := by
  ext σ
  simp [AddSubgroup.mem_inertia]

open IsLocalRing in
/-- `G_0` は極大イデアルの inertia、すなわち**古典的な惰性群**(剰余体に自明に作用する
自己同型のなす部分群)。★`lowerRamificationGroup` の定義が `𝔪^{n+1}` なので `n = 0` で
`𝔪^1 = 𝔪`。 -/
theorem lowerRamificationGroup_zero_eq_inertia (B : Type*) [CommRing B] [IsLocalRing B]
    (G : Type*) [Group G] [MulSemiringAction G B] :
    lowerRamificationGroup B G 0 = (maximalIdeal B).inertia G := by
  rw [lowerRamificationGroup, zero_add, pow_one]

/-- `G_n ≤ G_0`(`Antitone` の `0 ≤ n` での値)。 -/
theorem lowerRamificationGroup_le_zero {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (n : ℕ) :
    lowerRamificationGroup B G n ≤ lowerRamificationGroup B G 0 :=
  lowerRamificationGroup_antitone B G (Nat.zero_le n)

open IsLocalRing in
/-- 惰性群を含む部分群は、すべての `G_n`(`n ≥ 0`)を含む。 -/
theorem lowerRamificationGroup_le_of_inertia_le {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G}
    (h : (maximalIdeal B).inertia G ≤ H) (n : ℕ) :
    lowerRamificationGroup B G n ≤ H := by
  refine le_trans (lowerRamificationGroup_le_zero n) ?_
  rw [lowerRamificationGroup, zero_add, pow_one]
  exact h

/-- ★**抽象核(§1)**: `G_n ≤ H` なら、`H` の側で作った `G_n` は
`H.subtype` で `G` の側の `G_n` の上に写る。

★`(G_n).subgroupOf H = G_n(H)` が `rfl`(`subgroupOf_lowerRamificationGroup`)なので、
中身は mathlib の `Subgroup.map_subgroupOf_eq_of_le` だけである。 -/
theorem map_subtype_lowerRamificationGroup {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} {n : ℕ}
    (h : lowerRamificationGroup B G n ≤ H) :
    Subgroup.map H.subtype (lowerRamificationGroup B (H : Type _) n)
      = lowerRamificationGroup B G n := by
  rw [← subgroupOf_lowerRamificationGroup, Subgroup.map_subgroupOf_eq_of_le h]

open IsLocalRing in
/-- ★★★**§1 の本体 —— 不分岐底変換で下付き分岐群は変わらない**。

`B = 𝒪_L` を固定し、群だけを `G = Gal(L/K)` から `H = Gal(L/L₀)` へ縮める。
惰性群 `𝔪.inertia G`(`= G_0`)が `H` に入っていれば、すべての `n ≥ 0` で
`G_n(L/L₀)` は `G_n(L/K)` に一致する。

★「付値が同じという理由だけ」の正体がこれである——`L₀` は `L` の中の話ではなく
**底**の話なので、整数環も極大イデアルもまったく変わらない。
変わるのは作用する群だけで、`G_n ⊆ G_0 ⊆ H` だから群の縮小が `G_n` に効かない。 -/
theorem map_subtype_lowerRamificationGroup_of_inertia_le
    {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G}
    (h : (maximalIdeal B).inertia G ≤ H) (n : ℕ) :
    Subgroup.map H.subtype (lowerRamificationGroup B (H : Type _) n)
      = lowerRamificationGroup B G n :=
  map_subtype_lowerRamificationGroup (lowerRamificationGroup_le_of_inertia_le h n)

/-- 元の言葉での形。★`rfl`(`lowerRamificationGroup` が `Ideal.inertia` で定義されており、
`H` の元への制限が定義通りだから)。 -/
theorem mem_lowerRamificationGroup_subtype_iff {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} {n : ℕ} {σ : H} :
    σ ∈ lowerRamificationGroup B (H : Type _) n ↔ (σ : G) ∈ lowerRamificationGroup B G n :=
  Iff.rfl

/-- 群同型としての形 `G_n(H) ≃* G_n(G)`。 -/
noncomputable def lowerRamificationGroupMulEquiv {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} {n : ℕ}
    (h : lowerRamificationGroup B G n ≤ H) :
    lowerRamificationGroup B (H : Type _) n ≃* lowerRamificationGroup B G n :=
  (Subgroup.equivMapOfInjective _ H.subtype Subtype.val_injective).trans
    (MulEquiv.subgroupCongr (map_subtype_lowerRamificationGroup h))

/-! ## §2 抽象核 —— 固定部分群の像は射影系と両立する

★★純群論。ここも分岐・付値・Galois の語彙が 1 つも出てこない。

Y19 の `StageFiltration.compat` は「`N ≤ M` のとき `S N v · M = S M v`」という形をしている。
`S N v := I ⊔ N`(`I` は `v` に依らない固定部分群)がこれを満たす、というのが本節である。 -/

/-- ★**核に入る部分を足しても像は変わらない**。 -/
theorem map_sup_of_le_ker {Γ G : Type*} [Group Γ] [Group G] (f : Γ →* G)
    (I N' : Subgroup Γ) (h : N' ≤ f.ker) : Subgroup.map f (I ⊔ N') = Subgroup.map f I := by
  rw [Subgroup.map_sup, sup_eq_left, Subgroup.map_le_iff_le_comap, Subgroup.comap_map_eq]
  exact le_trans h le_sup_right

/-- ★★**`StageFiltration.compat` の形**: `N′ ≤ N`(`N` は正規)なら
`(I ⊔ N′) · N = I ⊔ N`。 -/
theorem coe_sup_mul_coe_eq {Γ : Type*} [Group Γ] (I N' N : Subgroup Γ) [N.Normal] (h : N' ≤ N) :
    ((I ⊔ N' : Subgroup Γ) : Set Γ) * (N : Set Γ) = ((I ⊔ N : Subgroup Γ) : Set Γ) := by
  rw [← Subgroup.mul_normal, sup_assoc, sup_eq_right.2 h]

/-- `f` を `I ⊔ N′` に制限して像 `f(I)` へ落とした準同型(`N′ ≤ ker f`)。 -/
def restrictToMap {Γ G : Type*} [Group Γ] [Group G] (f : Γ →* G) (I N' : Subgroup Γ)
    (h : N' ≤ f.ker) : (↥(I ⊔ N') : Type _) →* (↥(Subgroup.map f I) : Type _) :=
  (f.comp (I ⊔ N').subtype).codRestrict _ (fun x => by
    have hx : f x.1 ∈ Subgroup.map f (I ⊔ N') := ⟨x.1, x.2, rfl⟩
    rwa [map_sup_of_le_ker f I N' h] at hx)

/-- ★**全射性**(§2 の中身)。 -/
theorem surjective_restrictToMap {Γ G : Type*} [Group Γ] [Group G] (f : Γ →* G) (I N' : Subgroup Γ)
    (h : N' ≤ f.ker) : Function.Surjective (restrictToMap f I N' h) := by
  rintro ⟨y, hy⟩
  obtain ⟨x, hx, rfl⟩ := hy
  exact ⟨⟨x, Subgroup.mem_sup_left hx⟩, rfl⟩

/-- ★**その核**は `ker f`(を `I ⊔ N′` の中で見たもの)。 -/
theorem ker_restrictToMap {Γ G : Type*} [Group Γ] [Group G] (f : Γ →* G) (I N' : Subgroup Γ)
    (h : N' ≤ f.ker) : (restrictToMap f I N' h).ker = f.ker.subgroupOf (I ⊔ N') := by
  ext x
  simp [restrictToMap, MonoidHom.mem_ker, Subgroup.mem_subgroupOf]

/-- ★★**§2 抽象核の到達点**: `(I ⊔ N′) / ker f ≅ f(I)`。

具体層では `Γ = Γ_K`、`f = (Γ_K ↠ Gal(L/K))`、`I = I_K`、`N′ = Gal(K̄/L′)` であり、
左辺が `Gal(L′/L′₀)`、右辺が `Gal(L/L₀)` に対応する。 -/
noncomputable def quotientEquivMap {Γ G : Type*} [Group Γ] [Group G] (f : Γ →* G)
    (I N' : Subgroup Γ) (h : N' ≤ f.ker) :
    (↥(I ⊔ N') ⧸ f.ker.subgroupOf (I ⊔ N')) ≃* ↥(Subgroup.map f I) :=
  (QuotientGroup.quotientMulEquivOfEq (ker_restrictToMap f I N' h).symm).trans
    (QuotientGroup.quotientKerEquivOfSurjective _ (surjective_restrictToMap f I N' h))

/-! ## §3 具体層(§1 の橋) —— `G_0` は剰余体への還元の核

`Ideal.inertia` で定義した `G_0` が、古典的な惰性群
`ker(Gal(K(x)/K) → Gal(k_{K(x)}/k))` に一致することの確認。 -/

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NNReal Valued

variable {p : ℕ} [Fact p.Prime]

/-- ★★**`G_0 = ker(還元射)`**——`Ideal.inertia` による定義が古典的な惰性群と一致する。

★これで §1 の仮定 `𝔪.inertia G ≤ H` は
「`H` は剰余体に自明に作用する自己同型をすべて含む」と読める。 -/
theorem lowerRamificationGroupAdjoin_zero_eq_ker_residueGalHom (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    lowerRamificationGroupAdjoin K x 0 = (residueGalHom K x).ker := by
  ext σ
  rw [MonoidHom.mem_ker, lowerRamificationGroupAdjoin, mem_lowerRamificationGroup_iff_forall]
  simp only [zero_add, pow_one]
  constructor
  · intro hσ
    refine AlgEquiv.ext fun z => ?_
    obtain ⟨b, rfl⟩ := IsLocalRing.residue_surjective z
    show residueAlgEquiv K x σ _ = _
    rw [residueAlgEquiv_apply, AlgEquiv.one_apply, ← sub_eq_zero, ← map_sub,
      IsLocalRing.residue_eq_zero_iff]
    exact hσ b
  · intro hσ b
    rw [← IsLocalRing.residue_eq_zero_iff, map_sub, sub_eq_zero]
    have h1 : residueAlgEquiv K x σ (IsLocalRing.residue (adjoinIntegers K x) b)
        = IsLocalRing.residue (adjoinIntegers K x) b := by
      show residueGalHom K x σ _ = _
      rw [hσ]; rfl
    rw [residueAlgEquiv_apply] at h1
    exact h1

/-! ## §4 具体層(§2) —— 惰性部分群 `Gal(L/L₀)` とその全射性 -/

/-- **絶対惰性群** `I_K = Gal(K̄/K^ur)`。

★`Found/PGC/InertiaIdentification.lean` の `inertia_eq_fixingSubgroup_unramifiedClosure` に
より、これは `Skeleton/PGC/Section1Defs.lean` が構成した `inertia` と同じものである
(`absInertia_eq_inertia`)。 -/
noncomputable def absInertia (K : PAdicLocalField p) : Subgroup K.absGal :=
  (unramifiedClosure K).fixingSubgroup

theorem absInertia_eq_inertia (K : PAdicLocalField p) :
    absInertia K = inertia (residueCardinality p) (subgroupCorrespondence p) K :=
  (inertia_eq_fixingSubgroup_unramifiedClosure K).symm

/-- `I_K` は正規部分群(`K^ur/K` が Galois だから)。 -/
instance normal_absInertia (K : PAdicLocalField p) : (absInertia K).Normal := by
  haveI := isGalois_closure K
  exact (InfiniteGalois.normal_iff_isGalois _).mpr (isGalois_unramifiedClosure K)

/-- `I_K` は閉部分群(固定化部分群は常に閉)。 -/
theorem isClosed_absInertia (K : PAdicLocalField p) :
    IsClosed ((absInertia K : Subgroup K.absGal) : Set K.absGal) :=
  InfiniteGalois.fixingSubgroup_isClosed (unramifiedClosure K)

/-- **有限段の惰性部分群** `I(L/K) ⊆ Gal(L/K)` ——絶対惰性群 `I_K` の像。

★★これが `Gal(L/L₀)`(`L₀ = L ⊓ K^ur`)であることは
`lift_fixedField_inertiaGal` + `IntermediateField.fixingSubgroup_fixedField` による。 -/
noncomputable def inertiaGal (K : PAdicLocalField p) (L : IntermediateField K.carrier K.closure)
    [Normal K.carrier L] : Subgroup (L ≃ₐ[K.carrier] L) :=
  Subgroup.map (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _))
    (absInertia K)

/-- ★★★**`L₀ = L ⊓ K^ur` は `inertiaGal K L` の固定体**。

mathlib の `InfiniteGalois.restrict_fixedField`(`fixedField H ⊓ L = lift (fixedField (H の像))`)
に `H := I_K` を代入し、`fixedField I_K = K^ur` を使うだけ。

★★これが「`L₀ := L ∩ K^ur`」と「`Gal(L/L₀) := I_K の像`」が同じものだ、という
本持ち場の中心の等式である。 -/
theorem lift_fixedField_inertiaGal (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    IntermediateField.lift (IntermediateField.fixedField (inertiaGal K L))
      = unramifiedClosure K ⊓ L := by
  haveI := isGalois_closure K
  rw [inertiaGal, absInertia,
    ← InfiniteGalois.restrict_fixedField ((unramifiedClosure K).fixingSubgroup) L,
    InfiniteGalois.fixedField_fixingSubgroup]

/-- `inertiaGal K L` は自分の固定体の固定化部分群(有限次の Galois 対応)。
★`lift_fixedField_inertiaGal` と合わせると `inertiaGal K L = Gal(L/L₀)`。 -/
theorem inertiaGal_eq_fixingSubgroup (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    [FiniteDimensional K.carrier L] :
    inertiaGal K L = (IntermediateField.fixedField (inertiaGal K L)).fixingSubgroup :=
  (IntermediateField.fixingSubgroup_fixedField _).symm

/-- `L/K` が normal なら `Gal(K̄/L)` は `Γ_K` の正規部分群(制限射の核だから)。
★`Subgroup.Normal` のインスタンスとして置いておかないと、§2 の商群を**式に書けない**
(`lean-idioms.md` の「期待型が無い `haveI` は遅い」の裏返しで、
型に現れるインスタンスは証明の中で `haveI` しても間に合わない)。 -/
instance normal_fixingSubgroup_of_normal (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    (L.fixingSubgroup : Subgroup K.absGal).Normal := by
  rw [← IntermediateField.restrictNormalHom_ker
    (K := K.carrier) (L := K.closure) (E := L)]
  infer_instance

/-- ★**`Γ_K` への引き戻し**: `Gal(L/L₀)` の逆像は `I_K ⊔ Gal(K̄/L)`。

★これで §2 の抽象核(`I ⊔ N` の形)がそのまま使える。 -/
theorem comap_inertiaGal (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    Subgroup.comap (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _))
        (inertiaGal K L)
      = absInertia K ⊔ L.fixingSubgroup := by
  rw [inertiaGal, Subgroup.comap_map_eq, IntermediateField.restrictNormalHom_ker]

/-- ★★★**§2 の本体(全射性)**: `L ⊆ L′` のとき
`Gal(L′/L′₀)`(の `Γ_K` での姿 `I_K ⊔ Gal(K̄/L′)`)は `Gal(L/L₀)` の上へ写る。

★中身は「核に入る部分を足しても像は変わらない」(`map_sup_of_le_ker`)だけである
——`Gal(K̄/L′) ≤ Gal(K̄/L) = ker(Γ_K ↠ Gal(L/K))`。
★★**同じ `I_K` の像**を取っているので、全射性は自動的に出る。
Y19 が測った「道 A/B」が正しかったことの確認でもある。 -/
theorem map_restrictNormalHom_sup_fixingSubgroup (K : PAdicLocalField p)
    {L L' : IntermediateField K.carrier K.closure} [Normal K.carrier L] (hLL' : L ≤ L') :
    Subgroup.map (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _))
        (absInertia K ⊔ L'.fixingSubgroup)
      = inertiaGal K L := by
  rw [inertiaGal]
  refine map_sup_of_le_ker _ _ _ ?_
  rw [IntermediateField.restrictNormalHom_ker]
  exact IntermediateField.fixingSubgroup_le hLL'

/-- ★**その核**: `Gal(K̄/L)` を `I_K ⊔ Gal(K̄/L′)` の中で見たもの。

有限段に降ろすと `Gal(L′/L′₀) ⊓ Gal(L′/L) = Gal(L′/L·L′₀)` であり、
Corollary 6.13(i) の `H` に対応する(ファイル冒頭「退化の自己検査」参照)。 -/
noncomputable def inertiaQuotientEquiv (K : PAdicLocalField p)
    {L L' : IntermediateField K.carrier K.closure} [Normal K.carrier L] (hLL' : L ≤ L') :
    (↥(absInertia K ⊔ L'.fixingSubgroup) ⧸
        L.fixingSubgroup.subgroupOf (absInertia K ⊔ L'.fixingSubgroup))
      ≃* ↥(inertiaGal K L) :=
  (QuotientGroup.quotientMulEquivOfEq (by
    rw [IntermediateField.restrictNormalHom_ker
      (K := K.carrier) (L := K.closure) (E := L)])).trans
    ((quotientEquivMap (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        (L : Type _)) (absInertia K) L'.fixingSubgroup (by
        rw [IntermediateField.restrictNormalHom_ker]
        exact IntermediateField.fixingSubgroup_le hLL')).trans
      (MulEquiv.subgroupCongr rfl))

/-- ★★**`compat`(Y19 の `StageFiltration` が要求する形)**:
`L ⊆ L′` のとき `(I_K ⊔ Gal(K̄/L′)) · Gal(K̄/L) = I_K ⊔ Gal(K̄/L)`。 -/
theorem coe_sup_fixingSubgroup_mul_coe (K : PAdicLocalField p)
    {L L' : IntermediateField K.carrier K.closure} [Normal K.carrier L] (hLL' : L ≤ L') :
    ((absInertia K ⊔ L'.fixingSubgroup : Subgroup K.absGal) : Set K.absGal)
        * ((L.fixingSubgroup : Subgroup K.absGal) : Set K.absGal)
      = ((absInertia K ⊔ L.fixingSubgroup : Subgroup K.absGal) : Set K.absGal) := by
  exact coe_sup_mul_coe_eq _ _ _ (IntermediateField.fixingSubgroup_le hLL')

/-! ## §5 整合性検査 —— Yoshida §6.1 の設定(完全分岐)では `G_0 = Gal(L/L₀)`

★一般形は**埋まっていない**(ファイル冒頭「入っていないもの」)。ここで確かめるのは
原典が実際に扱っている場合、すなわち `L/K` が完全分岐のときである。 -/

/-- `L ⊓ K^ur = ⊥` なら `inertiaGal K L = ⊤`(= `Gal(L/L₀)` で `L₀ = K`)。 -/
theorem inertiaGal_eq_top_of_inf_eq_bot (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    [FiniteDimensional K.carrier L] (h : unramifiedClosure K ⊓ L = ⊥) :
    inertiaGal K L = ⊤ := by
  have h1 : IntermediateField.fixedField (inertiaGal K L) = ⊥ := by
    rw [← IntermediateField.lift_inj, IntermediateField.lift_bot, lift_fixedField_inertiaGal, h]
  rw [inertiaGal_eq_fixingSubgroup K L, h1, IntermediateField.fixingSubgroup_bot]

/-- ★★**整合性検査**: `K(x)/K` が完全分岐なら `inertiaGal K K(x) = G_0`(どちらも `⊤`)。

★本体の見当「`G_0 = 惰性群 = Gal(L/L₀)`」が、原典が実際に扱っている
完全分岐の場合には**成り立っている**ことの確認。
★一般の有限 Galois `L/K` については未証明である(冒頭参照)。 -/
theorem inertiaGal_eq_lowerRamificationGroupAdjoin_zero_of_isTotallyRamified
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    inertiaGal K (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = lowerRamificationGroupAdjoin K x 0 := by
  rw [lowerRamificationGroupAdjoin_zero_eq_top K x ht,
    inertiaGal_eq_top_of_inf_eq_bot K _
      (by rw [inf_comm]; exact totallyRamifiedAdjoin_inf_unramifiedClosure K ht)]

/-- ★★★★**`Gal(L/L₀) ≤ G_0`**(= `inertiaGal K L ≤ G_0(L/K)`)。

★★これが本ファイルで一番中身のある定理である。証明は **Teichmüller 表現**:
`b ∈ 𝒪_L` の剰余が `0` でなければ、`ζ^{Q−1} = 1` かつ `ζ ≡ b` なる `ζ ∈ 𝒪_L` が取れる
(`ABC3.Found.exists_teichmullerRep`)。`p ∤ Q − 1` なので `ζ ∈ K^ur`
(`mem_unramifiedClosure_of_pow_eq_one`)であり、`σ` を延長した `τ ∈ I_K` は `ζ` を固定する。
よって `σ b − b = σ(b − ζ) + (ζ − b) ∈ 𝔪`。剰余が `0` の場合は
`smul_mem_maximalIdeal` だけで済む。

★**方向に注意**: これは `inertiaGal ≤ G_0` であって逆ではない。
逆(`G_0 ≤ inertiaGal`)は**埋まっていない**(ファイル冒頭)。 -/
theorem inertiaGal_le_lowerRamificationGroupAdjoin_zero (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    inertiaGal K (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≤ lowerRamificationGroupAdjoin K x 0 := by
  rintro σ ⟨τ, hτ, rfl⟩
  rw [absInertia] at hτ
  have hτ0 : τ ∈ (unramifiedClosure K).fixingSubgroup := hτ
  rw [IntermediateField.mem_fixingSubgroup_iff] at hτ0
  rw [lowerRamificationGroupAdjoin, mem_lowerRamificationGroup_iff_forall]
  simp only [zero_add, pow_one]
  intro b
  set L := IntermediateField.adjoin K.carrier ({x} : Set K.closure) with hL
  set σ := AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (L : Type _) τ with hσdef
  by_cases hb : IsLocalRing.residue (adjoinIntegers K x) b = 0
  · have hbm : b ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) := by
      rw [← IsLocalRing.residue_eq_zero_iff]; exact hb
    exact Ideal.sub_mem _ (smul_mem_maximalIdeal σ hbm) hbm
  · haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K x)) := Fintype.ofFinite _
    set Q := Fintype.card (IsLocalRing.ResidueField (adjoinIntegers K x)) with hQdef
    have hQ2 : 2 ≤ Q := Fintype.one_lt_card
    have hpQ : p ∣ Q := by
      rw [hQdef, ← Nat.card_eq_fintype_card]
      exact prime_dvd_card_residueField_adjoinIntegers K x
    have hpQ1 : ¬ p ∣ (Q - 1) := by
      intro hcon
      have h1 : p ∣ (Q - (Q - 1)) := Nat.dvd_sub hpQ hcon
      rw [show Q - (Q - 1) = 1 by omega] at h1
      exact (Nat.Prime.one_lt (Fact.out : p.Prime)).ne' (Nat.dvd_one.mp h1)
    obtain ⟨ζ, hζpow, hζres⟩ :=
      ABC3.Found.exists_teichmullerRep (IsLocalRing.residue (adjoinIntegers K x) b) hb
    set ι : adjoinIntegers K x →+* K.closure :=
      (algebraMap (L : Type _) K.closure).comp (SubringClass.subtype (adjoinIntegers K x))
      with hιdef
    have hιinj : Function.Injective ι := by
      intro u v huv
      exact Subtype.ext ((algebraMap (L : Type _) K.closure).injective huv)
    have hmem : ι ζ ∈ unramifiedClosure K :=
      mem_unramifiedClosure_of_pow_eq_one K (by omega) hpQ1
        (by rw [← map_pow, hζpow, map_one])
    have hfix : τ (ι ζ) = ι ζ := hτ0 _ hmem
    have hσζ : σ • ζ = ζ := by
      apply hιinj
      show algebraMap (L : Type _) K.closure ((σ • ζ : adjoinIntegers K x) : (L : Type _))
        = ι ζ
      rw [coe_smul_adjoinIntegers]
      have hres := AlgEquiv.restrictNormalHom_apply (F := K.carrier) (K₁ := K.closure)
        L τ ((ζ : adjoinIntegers K x) : (L : Type _))
      show ((σ ((ζ : adjoinIntegers K x) : (L : Type _)) : (L : Type _)) : K.closure) = _
      rw [hσdef, hres]
      exact hfix
    have hbz : b - ζ ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) := by
      rw [← IsLocalRing.residue_eq_zero_iff, map_sub, hζres, sub_self]
    have hsplit : σ • b - b = σ • (b - ζ) + (ζ - b) := by
      rw [smul_sub, hσζ]; ring
    rw [hsplit]
    exact Ideal.add_mem _ (smul_mem_maximalIdeal σ hbz)
      (by rw [← neg_sub]; exact neg_mem hbz)

/-- ★★★★**`L/L₀` は完全分岐**(Yoshida §6.1 が要求する `G = G_0` の形)。

`G_0(L/L₀) = (G_0(L/K)).subgroupOf Gal(L/L₀) = ⊤`。★これが本持ち場の
「不分岐部分の除去」の到達点である——一般の有限 Galois `L/K` を、
不分岐部分 `L₀ = L ⊓ K^ur` へ底を移すことで**原典の設定に落とせる**。

★★ここで §1 の抽象核が効いている: `(G_0).subgroupOf H = G_0(H)` は `rfl`
(`subgroupOf_lowerRamificationGroupAdjoin`)なので、必要なのは包含
`H ≤ G_0` だけである。 -/
theorem lowerRamificationGroup_inertiaGal_zero_eq_top (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    lowerRamificationGroup (adjoinIntegers K x)
        (inertiaGal K (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : Type _) 0
      = ⊤ := by
  rw [← subgroupOf_lowerRamificationGroupAdjoin, Subgroup.subgroupOf_eq_top]
  exact inertiaGal_le_lowerRamificationGroupAdjoin_zero K x

open IsLocalRing in
/-- ★**§1 を具体層へ**: `Gal(K(x)/K)` の部分群 `H` が惰性群 `G_0` を含めば、
`H` の側で作った `G_n` は `G_n` そのものである。

★仮定 `hH` は「`H ⊇ Gal(K(x)/K(x) ⊓ K^ur)`」であってほしいが、
その同定(`G_0 = inertiaGal`)は本ファイルの外である(冒頭「入っていないもの」)。
ここでは仮定のまま受け取る。 -/
theorem lowerRamificationGroupAdjoin_map_subtype (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {H : Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))}
    (hH : lowerRamificationGroupAdjoin K x 0 ≤ H) (n : ℕ) :
    Subgroup.map H.subtype (lowerRamificationGroup (adjoinIntegers K x) (H : Type _) n)
      = lowerRamificationGroupAdjoin K x n :=
  map_subtype_lowerRamificationGroup (le_trans (lowerRamificationGroup_le_zero n) hH)

/-! ## §6 `I_K` から作った段データ(★**原典の `Γ_K^v` ではない**)

★★**警告**: 以下の `inertiaStageFiltration` は `v ≥ 0` で `Γ^v = I_K` を与える。
真の `Γ_K^v` は `v > 0` で野性惰性群に含まれる(もっと小さい)ので、**これは上からの近似**
である。Y19 の `trivialStageFiltration`(`v > 0` で `⊥`)が下からの退化 witness だったのに
対応する。★`Skeleton/PGC/Section2.lean` の入力にも G2 の非空虚 witness にも
**使ってはならない**。

置いてある理由は 1 つだけ: §2 の `compat` が Y19 の `StageFiltration.ofOpenNormal` に
**実際に嵌まる**ことを型で検査するためである。 -/

/-- `I_K` を使った段データ。`v < 0` で `⊤`、`v ≥ 0` で `I_K ⊔ N`。
★**原典の `Γ_K^v` ではない**(節冒頭の警告)。 -/
noncomputable def inertiaStageFiltration (K : PAdicLocalField p) : StageFiltration K.absGal :=
  StageFiltration.ofOpenNormal (fun N v => if v < 0 then ⊤ else absInertia K ⊔ N)
    (fun _ _ v => by
      by_cases hv : v < 0
      · simp [hv]
      · simp only [hv, if_false]; exact le_sup_right)
    (fun _ hN v => by
      haveI := hN.2
      by_cases hv : v < 0
      · simp only [hv, if_true]; infer_instance
      · simp only [hv, if_false]; infer_instance)
    (fun _ _ a b hab => by
      by_cases hb : b < 0
      · have ha : a < 0 := lt_of_le_of_lt hab hb
        simp [ha, hb]
      · by_cases ha : a < 0 <;> simp [ha, hb])
    (fun M hM _ _ hNM v => by
      haveI := hM.2
      by_cases hv : v < 0
      · simp only [hv, if_true]
        rw [← Subgroup.mul_normal, top_sup_eq]
      · simp only [hv, if_false]
        exact coe_sup_mul_coe_eq _ _ _ hNM)

/-- `v ≥ 0` での極限は `I_K`。★`I_K` が閉であることと、開正規部分群が `1` の近傍基を
なすことから、Y19 の一意性定理 `eq_limit_of_isClosed` で出る。 -/
theorem inertiaStageFiltration_limit_of_nonneg (K : PAdicLocalField p) {v : ℝ} (hv : ¬ v < 0) :
    (inertiaStageFiltration K).limit v = absInertia K := by
  refine ((inertiaStageFiltration K).eq_limit_of_isClosed ?_ (isClosed_absInertia K) ?_).symm
  · intro U hU
    obtain ⟨N, hN, hNU⟩ := absGal_exists_mem_openNormalBase_subset K hU
    exact ⟨N, hN, hNU⟩
  · intro M hM
    haveI := hM.2
    show _ = ((if v < 0 then (⊤ : Subgroup K.absGal) else absInertia K ⊔ M : Subgroup K.absGal)
      : Set K.absGal)
    simp only [hv, if_false]
    rw [← Subgroup.mul_normal]

/-- `v < 0` での極限は `⊤`。 -/
theorem inertiaStageFiltration_limit_of_neg (K : PAdicLocalField p) {v : ℝ} (hv : v < 0) :
    (inertiaStageFiltration K).limit v = ⊤ := by
  rw [eq_top_iff]
  intro x _
  rw [StageFiltration.mem_limit_iff]
  intro N _
  simp [inertiaStageFiltration, hv]

end ABC3.Found.PGC
